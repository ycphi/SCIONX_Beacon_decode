"""Multi-phase soft HDLC decoder with CRC-syndrome repair.

The decoder deliberately separates four decisions:

1. demodulate the same burst at every samples-per-symbol phase;
2. use a bounded beam search to recover a legal HDLC stuffing path and the
   configured minimum number of closing flags (three by default);
3. repair the first non-empty CRC syndrome at Hamming weight 1..4; and
4. accept a result only when all phases/paths collapse to one distinct frame.

The hard 16-byte SCION-X AX.25 address/control/PID prefix and the exact frame
length are protocol constraints.  CRC is the final full-frame validator.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from typing import Sequence

import numpy as np

from .auto_repair import (
    EXPECTED_FRAME_BITS,
    FLAG_BITS,
    HARD_ANCHOR,
    PAYLOAD_BYTES,
    bits_lsb_first,
    bits_to_bytes,
    stuff_bits,
    strict_destuff,
)
from .hdlc import crc16_x25


@dataclass(frozen=True)
class PhaseObservation:
    phase: int
    offset: int
    raw_bits: np.ndarray
    soft_distance: np.ndarray
    header_errors: int = 0
    verifiable_mismatches: int = 0
    verifiable_count: int = 0


@dataclass(frozen=True)
class _BeamState:
    output_count: int
    ones_run: int
    flip_indices: tuple[int, ...]
    cost: float


@dataclass(frozen=True)
class _StructuralCandidate:
    phase: int
    offset: int
    frame: bytes
    flag_start: int
    flip_indices: tuple[int, ...]
    excluded_indices: tuple[int, ...]
    source_indices: tuple[int, ...]
    cost: float


@dataclass(frozen=True)
class _ValidatedSolution:
    candidate: _StructuralCandidate
    syndrome_destuffed_indices: tuple[int, ...]
    syndrome_raw_indices: tuple[int, ...]
    frame: bytes
    cost: float

    @property
    def all_raw_flips(self) -> tuple[int, ...]:
        return tuple(sorted(set(
            self.candidate.flip_indices + self.syndrome_raw_indices
        )))


@dataclass(frozen=True)
class AdvancedRepairResult:
    status: str
    reason: str
    original_stuffing_aligned: bool
    original_crc_pass: bool
    hard_anchor_match: bool
    closing_flag_found: bool
    selected_flag_start: int | None
    destuffed_bit_count: int
    excluded_indices: tuple[int, ...]
    flip_indices: tuple[int, ...]
    source_indices: tuple[int, ...]
    candidate_pool: tuple[int, ...]
    candidate_evaluations: int
    crc_solution_count: int
    repair_cost: float | None
    frame: bytes | None
    selected_phase: int | None
    selected_offset: int | None
    phases_tested: tuple[int, ...]
    beam_candidate_count: int
    beam_flip_indices: tuple[int, ...]
    syndrome_flip_destuffed_indices: tuple[int, ...]
    syndrome_flip_raw_indices: tuple[int, ...]
    syndrome_flip_count: int
    arbitration_distinct_frames: int
    arbitration_score: float | None
    arbitration_margin: float | None
    repair_method: str

    @property
    def auto_repaired(self) -> bool:
        return self.status == "AUTO_REPAIRED_UNIQUE"

    def to_summary(self) -> dict:
        return {
            "status": self.status,
            "reason": self.reason,
            "original_stuffing_aligned": self.original_stuffing_aligned,
            "original_crc_pass": self.original_crc_pass,
            "hard_anchor_match": self.hard_anchor_match,
            "closing_flag_found": self.closing_flag_found,
            "selected_flag_start": self.selected_flag_start,
            "destuffed_bit_count": self.destuffed_bit_count,
            "excluded_indices": list(self.excluded_indices),
            "flip_indices": list(self.flip_indices),
            "candidate_pool": list(self.candidate_pool),
            "candidate_evaluations": self.candidate_evaluations,
            "crc_solution_count": self.crc_solution_count,
            "repair_cost": self.repair_cost,
            "frame_hex": self.frame.hex() if self.frame is not None else "",
            "selected_phase": self.selected_phase,
            "selected_offset": self.selected_offset,
            "phases_tested": list(self.phases_tested),
            "beam_candidate_count": self.beam_candidate_count,
            "beam_flip_indices": list(self.beam_flip_indices),
            "syndrome_flip_destuffed_indices": list(
                self.syndrome_flip_destuffed_indices
            ),
            "syndrome_flip_raw_indices": list(self.syndrome_flip_raw_indices),
            "syndrome_flip_count": self.syndrome_flip_count,
            "arbitration_distinct_frames": self.arbitration_distinct_frames,
            "arbitration_score": self.arbitration_score,
            "arbitration_margin": self.arbitration_margin,
            "repair_method": self.repair_method,
        }


def _normalise_observation(
    observation: PhaseObservation,
    minimum_closing_flags: int,
) -> PhaseObservation:
    bits = np.asarray(observation.raw_bits, dtype=np.uint8).reshape(-1)
    distances = np.asarray(observation.soft_distance, dtype=np.float64).reshape(-1)
    if bits.size != distances.size:
        raise ValueError("Each phase must have equal raw_bits and soft_distance lengths.")
    if bits.size < EXPECTED_FRAME_BITS + len(FLAG_BITS) * minimum_closing_flags:
        raise ValueError("A phase observation is too short for one SCION-X frame.")
    distances = np.nan_to_num(
        distances, nan=np.inf, posinf=np.inf, neginf=np.inf
    )
    return PhaseObservation(
        phase=int(observation.phase),
        offset=int(observation.offset),
        raw_bits=bits,
        soft_distance=distances,
        header_errors=int(observation.header_errors),
        verifiable_mismatches=int(observation.verifiable_mismatches),
        verifiable_count=int(observation.verifiable_count),
    )


def _candidate_pool(bits: np.ndarray, distances: np.ndarray, size: int) -> tuple[int, ...]:
    if size <= 0:
        return tuple()
    limit = min(bits.size, EXPECTED_FRAME_BITS + 160)
    anchor = bits_lsb_first(HARD_ANCHOR)
    mandatory = {
        index for index in range(min(len(anchor), limit))
        if int(bits[index]) != int(anchor[index])
    }
    protocol: set[int] = set()
    for index in range(5, limit):
        previous = bits[index - 5:index]
        ones = int(np.count_nonzero(previous))
        if ones == 5 and int(bits[index]) == 1:
            protocol.update(range(index - 5, index + 1))
        elif ones == 4 and int(bits[index]) == 0:
            protocol.add(index)
    ranked = sorted(
        range(limit),
        key=lambda index: (
            0 if index in mandatory else 1,
            float(distances[index]),
            0 if index in protocol else 1,
            index,
        ),
    )
    return tuple(sorted(ranked[:max(size, len(mandatory))]))


def _add_flip(
    state: _BeamState,
    raw_index: int,
    distance: float,
    maximum_flips: int,
) -> tuple[tuple[int, ...], float] | None:
    if raw_index in state.flip_indices:
        return state.flip_indices, state.cost
    if len(state.flip_indices) >= maximum_flips:
        return None
    return state.flip_indices + (raw_index,), state.cost + float(distance)


def _finish_candidate(
    observation: PhaseObservation,
    state: _BeamState,
    next_raw_index: int,
    maximum_flips: int,
    minimum_closing_flags: int,
) -> _StructuralCandidate | None:
    bits = observation.raw_bits
    distances = observation.soft_distance
    flips = state.flip_indices
    cost = state.cost
    raw_index = next_raw_index

    if state.ones_run == 5:
        if raw_index >= bits.size:
            return None
        if int(bits[raw_index]) != 0:
            added = _add_flip(
                _BeamState(state.output_count, state.ones_run, flips, cost),
                raw_index, distances[raw_index], maximum_flips,
            )
            if added is None:
                return None
            flips, cost = added
        raw_index += 1

    flag_start = raw_index
    expected_flags = FLAG_BITS * minimum_closing_flags
    if flag_start + len(expected_flags) > bits.size:
        return None
    for offset, expected in enumerate(expected_flags):
        index = flag_start + offset
        actual = int(bits[index]) ^ (index in flips)
        if actual == expected:
            continue
        temp_state = _BeamState(state.output_count, state.ones_run, flips, cost)
        added = _add_flip(temp_state, index, distances[index], maximum_flips)
        if added is None:
            return None
        flips, cost = added

    corrected = bits.copy()
    if flips:
        corrected[list(flips)] ^= 1
    destuffed = strict_destuff(corrected[:flag_start])
    if not destuffed.valid or len(destuffed.bits) != EXPECTED_FRAME_BITS:
        return None
    frame = bits_to_bytes(destuffed.bits)
    if frame[:len(HARD_ANCHOR)] != HARD_ANCHOR:
        return None
    return _StructuralCandidate(
        phase=observation.phase,
        offset=observation.offset,
        frame=frame,
        flag_start=flag_start,
        flip_indices=tuple(sorted(flips)),
        excluded_indices=destuffed.excluded_indices,
        source_indices=destuffed.source_indices,
        cost=cost,
    )


def _beam_structural_candidates(
    observation: PhaseObservation,
    *,
    beam_width: int,
    branch_pool_size: int,
    maximum_structural_flips: int,
    maximum_candidates: int,
    minimum_closing_flags: int,
) -> tuple[list[_StructuralCandidate], tuple[int, ...], int]:
    bits = observation.raw_bits
    distances = observation.soft_distance
    pool = _candidate_pool(bits, distances, branch_pool_size)
    pool_set = set(pool)
    anchor_bits = bits_lsb_first(HARD_ANCHOR)
    states = [_BeamState(0, 0, tuple(), 0.0)]
    candidates: dict[tuple[bytes, int], _StructuralCandidate] = {}
    evaluations = 0

    for raw_index in range(bits.size):
        next_states: list[_BeamState] = []
        for state in states:
            if state.output_count == EXPECTED_FRAME_BITS:
                candidate = _finish_candidate(
                    observation,
                    state,
                    raw_index,
                    maximum_structural_flips,
                    minimum_closing_flags,
                )
                evaluations += 1
                if candidate is not None:
                    key = (candidate.frame, candidate.flag_start)
                    previous = candidates.get(key)
                    if previous is None or candidate.cost < previous.cost:
                        candidates[key] = candidate
                continue

            hard = int(bits[raw_index])
            if state.ones_run == 5:
                effective_values = (0,)
            elif state.output_count < len(anchor_bits):
                effective_values = (int(anchor_bits[state.output_count]),)
            elif raw_index in pool_set:
                effective_values = (hard, 1 - hard)
            else:
                effective_values = (hard,)

            for effective in effective_values:
                flips = state.flip_indices
                cost = state.cost
                if effective != hard:
                    added = _add_flip(
                        state, raw_index, distances[raw_index],
                        maximum_structural_flips,
                    )
                    if added is None:
                        continue
                    flips, cost = added

                if state.ones_run == 5:
                    next_state = _BeamState(
                        state.output_count, 0, flips, cost
                    )
                else:
                    next_state = _BeamState(
                        state.output_count + 1,
                        state.ones_run + 1 if effective else 0,
                        flips,
                        cost,
                    )
                next_states.append(next_state)
                evaluations += 1

        if not next_states:
            break
        # Exact duplicate histories can arise after deterministic stuffing.
        unique: dict[tuple[int, int, tuple[int, ...]], _BeamState] = {}
        for state in next_states:
            key = (state.output_count, state.ones_run, state.flip_indices)
            previous = unique.get(key)
            if previous is None or state.cost < previous.cost:
                unique[key] = state
        states = sorted(
            unique.values(),
            key=lambda state: (state.cost, len(state.flip_indices), state.flip_indices),
        )[:beam_width]
        if len(candidates) >= maximum_candidates and all(
            state.output_count >= EXPECTED_FRAME_BITS for state in states
        ):
            break

    ranked = sorted(
        candidates.values(),
        key=lambda item: (item.cost, len(item.flip_indices), item.flag_start),
    )[:maximum_candidates]
    return ranked, pool, evaluations


@lru_cache(maxsize=1)
def _crc_syndrome_deltas() -> tuple[int, ...]:
    payload = bytearray(PAYLOAD_BYTES)
    baseline = crc16_x25(payload)
    deltas: list[int] = []
    for bit_index in range(PAYLOAD_BYTES * 8):
        byte_index, within_byte = divmod(bit_index, 8)
        payload[byte_index] ^= 1 << within_byte
        deltas.append(crc16_x25(payload) ^ baseline)
        payload[byte_index] ^= 1 << within_byte
    deltas.extend(1 << bit_index for bit_index in range(16))
    return tuple(deltas)


def _minimal_syndrome_combinations(
    positions: Sequence[int],
    syndrome: int,
    maximum_bits: int,
    maximum_solutions: int = 2,
) -> tuple[list[tuple[int, ...]], int]:
    deltas = _crc_syndrome_deltas()
    positions = tuple(int(position) for position in positions)
    evaluations = 0

    one_map: dict[int, list[int]] = {}
    for position in positions:
        one_map.setdefault(deltas[position], []).append(position)
        evaluations += 1
    if maximum_bits >= 1 and syndrome in one_map:
        answers = [(position,) for position in one_map[syndrome]]
        return answers[:maximum_solutions], evaluations

    pair_map: dict[int, list[tuple[int, int]]] = {}
    for left_index in range(len(positions)):
        left = positions[left_index]
        for right in positions[left_index + 1:]:
            value = deltas[left] ^ deltas[right]
            pair_map.setdefault(value, []).append((left, right))
            evaluations += 1
    if maximum_bits >= 2 and syndrome in pair_map:
        return pair_map[syndrome][:maximum_solutions], evaluations

    if maximum_bits >= 3:
        answers: set[tuple[int, ...]] = set()
        for position in positions:
            target = syndrome ^ deltas[position]
            for pair in pair_map.get(target, ()): 
                evaluations += 1
                if position in pair:
                    continue
                answers.add(tuple(sorted((position,) + pair)))
                if len(answers) >= maximum_solutions:
                    return sorted(answers), evaluations
        if answers:
            return sorted(answers), evaluations

    if maximum_bits >= 4:
        answers = set()
        for pair_value, left_pairs in pair_map.items():
            right_pairs = pair_map.get(syndrome ^ pair_value, ())
            if not right_pairs:
                continue
            for left in left_pairs:
                for right in right_pairs:
                    evaluations += 1
                    combined = left + right
                    if len(set(combined)) != 4:
                        continue
                    answers.add(tuple(sorted(combined)))
                    if len(answers) >= maximum_solutions:
                        return sorted(answers), evaluations
        if answers:
            return sorted(answers), evaluations
    return [], evaluations


def _syndrome_solutions(
    candidate: _StructuralCandidate,
    observation: PhaseObservation,
    *,
    maximum_syndrome_bits: int,
    syndrome_pool_size: int,
) -> tuple[list[_ValidatedSolution], tuple[int, ...], int]:
    received_fcs = int.from_bytes(
        candidate.frame[PAYLOAD_BYTES:PAYLOAD_BYTES + 2], "little"
    )
    syndrome = crc16_x25(candidate.frame[:PAYLOAD_BYTES]) ^ received_fcs
    if syndrome == 0:
        return [
            _ValidatedSolution(candidate, tuple(), tuple(), candidate.frame, candidate.cost)
        ], tuple(), 1

    ranked_positions = sorted(
        range(len(candidate.source_indices)),
        key=lambda position: (
            float(observation.soft_distance[candidate.source_indices[position]]),
            position,
        ),
    )
    # The fixed AX.25/SCION-X anchor is not eligible for CRC-only alteration.
    ranked_positions = [
        position for position in ranked_positions
        if position >= len(HARD_ANCHOR) * 8
        and candidate.source_indices[position] not in candidate.flip_indices
    ][:syndrome_pool_size]
    pool = tuple(sorted(ranked_positions))
    combinations_found, evaluations = _minimal_syndrome_combinations(
        pool, syndrome, maximum_syndrome_bits
    )

    solutions: list[_ValidatedSolution] = []
    for destuffed_indices in combinations_found:
        corrected = bytearray(candidate.frame)
        for bit_index in destuffed_indices:
            byte_index, within_byte = divmod(bit_index, 8)
            corrected[byte_index] ^= 1 << within_byte
        corrected_frame = bytes(corrected)
        if corrected_frame[:len(HARD_ANCHOR)] != HARD_ANCHOR:
            continue
        corrected_fcs = int.from_bytes(
            corrected_frame[PAYLOAD_BYTES:PAYLOAD_BYTES + 2], "little"
        )
        if crc16_x25(corrected_frame[:PAYLOAD_BYTES]) != corrected_fcs:
            continue
        raw_indices = tuple(candidate.source_indices[index] for index in destuffed_indices)
        added_cost = sum(
            float(observation.soft_distance[index]) for index in raw_indices
        )
        solutions.append(_ValidatedSolution(
            candidate=candidate,
            syndrome_destuffed_indices=tuple(destuffed_indices),
            syndrome_raw_indices=raw_indices,
            frame=corrected_frame,
            cost=candidate.cost + added_cost,
        ))
    return solutions, pool, evaluations


def _multiphase_consensus_cost(
    frame: bytes, observations: Sequence[PhaseObservation]
) -> float:
    expected = np.asarray(stuff_bits(bits_lsb_first(frame)), dtype=np.uint8)
    total = 0.0
    used_phases = 0
    for observation in observations:
        count = min(expected.size, observation.raw_bits.size)
        if count <= 0:
            continue
        mismatch = observation.raw_bits[:count] != expected[:count]
        total += float(np.sum(observation.soft_distance[:count][mismatch]))
        used_phases += 1
    return total / max(used_phases, 1)


def advanced_repair_frame(
    observations: Sequence[PhaseObservation],
    *,
    beam_width: int = 128,
    beam_branch_pool_size: int = 48,
    maximum_structural_flips: int = 4,
    maximum_beam_candidates_per_phase: int = 6,
    maximum_syndrome_bits: int = 4,
    syndrome_pool_size: int = 48,
    minimum_arbitration_margin: float = 0.05,
    minimum_closing_flags: int = 3,
) -> AdvancedRepairResult:
    """Decode all sample phases and accept only one distinct CRC-valid frame."""
    if not observations:
        raise ValueError("At least one phase observation is required.")
    if beam_width <= 0 or maximum_beam_candidates_per_phase <= 0:
        raise ValueError("Beam limits must be positive.")
    if not 0 <= maximum_structural_flips <= 8:
        raise ValueError("maximum_structural_flips must be between 0 and 8.")
    if not 0 <= maximum_syndrome_bits <= 4:
        raise ValueError("maximum_syndrome_bits must be between 0 and 4.")
    if minimum_arbitration_margin < 0:
        raise ValueError("minimum_arbitration_margin must be non-negative.")
    if not 1 <= minimum_closing_flags <= 16:
        raise ValueError("minimum_closing_flags must be between 1 and 16.")

    phase_items = [
        _normalise_observation(item, minimum_closing_flags)
        for item in observations
    ]
    phase_by_number = {item.phase: item for item in phase_items}
    all_candidates: list[_StructuralCandidate] = []
    raw_pools: dict[int, tuple[int, ...]] = {}
    evaluations = 0
    for observation in phase_items:
        candidates, pool, count = _beam_structural_candidates(
            observation,
            beam_width=beam_width,
            branch_pool_size=beam_branch_pool_size,
            maximum_structural_flips=maximum_structural_flips,
            maximum_candidates=maximum_beam_candidates_per_phase,
            minimum_closing_flags=minimum_closing_flags,
        )
        all_candidates.extend(candidates)
        raw_pools[observation.phase] = pool
        evaluations += count

    no_flip_candidates = [
        candidate for candidate in all_candidates if not candidate.flip_indices
    ]
    original_aligned = bool(no_flip_candidates)
    original_crc = any(
        crc16_x25(candidate.frame[:PAYLOAD_BYTES]) == int.from_bytes(
            candidate.frame[PAYLOAD_BYTES:PAYLOAD_BYTES + 2], "little"
        )
        for candidate in no_flip_candidates
    )

    solutions_by_frame: dict[bytes, list[_ValidatedSolution]] = {}
    syndrome_pools: dict[int, tuple[int, ...]] = {}
    for candidate in all_candidates:
        observation = phase_by_number[candidate.phase]
        solutions, pool, count = _syndrome_solutions(
            candidate,
            observation,
            maximum_syndrome_bits=maximum_syndrome_bits,
            syndrome_pool_size=syndrome_pool_size,
        )
        syndrome_pools[candidate.phase] = pool
        evaluations += count
        for solution in solutions:
            solutions_by_frame.setdefault(solution.frame, []).append(solution)

    distinct_count = len(solutions_by_frame)
    minimum_flip_count = min(
        (
            len(solution.all_raw_flips)
            for solutions in solutions_by_frame.values()
            for solution in solutions
        ),
        default=None,
    )
    competitive_by_frame = {
        frame: [
            solution for solution in solutions
            if len(solution.all_raw_flips) == minimum_flip_count
        ]
        for frame, solutions in solutions_by_frame.items()
        if any(
            len(solution.all_raw_flips) == minimum_flip_count
            for solution in solutions
        )
    }
    competitive_count = len(competitive_by_frame)
    consensus_scores = {
        frame: _multiphase_consensus_cost(frame, phase_items)
        for frame in competitive_by_frame
    }
    ranked_frames = sorted(
        consensus_scores,
        key=lambda frame: (consensus_scores[frame], frame),
    )
    arbitration_margin = (
        consensus_scores[ranked_frames[1]] - consensus_scores[ranked_frames[0]]
        if len(ranked_frames) > 1 else None
    )
    winning_frame = (
        ranked_frames[0]
        if len(ranked_frames) == 1
        or (
            arbitration_margin is not None
            and arbitration_margin >= minimum_arbitration_margin
        )
        else None
    )
    phases = tuple(sorted(phase_by_number))
    if winning_frame is not None:
        same_frame_solutions = competitive_by_frame[winning_frame]
        selected = min(
            same_frame_solutions,
            key=lambda solution: (
                len(solution.all_raw_flips), solution.cost,
                solution.candidate.phase, solution.candidate.offset,
            ),
        )
        candidate = selected.candidate
        all_flips = selected.all_raw_flips
        status = "ORIGINAL_CRC_PASS" if not all_flips else "AUTO_REPAIRED_UNIQUE"
        method_parts = ["MULTIPHASE_HDLC_BEAM"]
        if selected.syndrome_destuffed_indices:
            method_parts.append(
                f"CRC_SYNDROME_{len(selected.syndrome_destuffed_indices)}BIT"
            )
        if distinct_count == 1:
            reason = (
                "All phase/path candidates collapse to one original CRC-valid frame."
                if status == "ORIGINAL_CRC_PASS"
                else "All phase/path candidates collapse to one CRC-valid repaired frame."
            )
        elif competitive_count == 1:
            reason = (
                f"{distinct_count} CRC-valid frames were found, but exactly one "
                f"has the minimum correction weight of {minimum_flip_count} bit(s)."
            )
        else:
            reason = (
                f"{competitive_count} minimum-weight CRC-valid frames were found; "
                f"multi-phase soft consensus selected one with margin "
                f"{arbitration_margin:.6g}."
            )
        return AdvancedRepairResult(
            status=status,
            reason=reason,
            original_stuffing_aligned=original_aligned,
            original_crc_pass=original_crc,
            hard_anchor_match=True,
            closing_flag_found=True,
            selected_flag_start=candidate.flag_start,
            destuffed_bit_count=EXPECTED_FRAME_BITS,
            excluded_indices=candidate.excluded_indices,
            flip_indices=all_flips,
            source_indices=candidate.source_indices,
            candidate_pool=raw_pools.get(candidate.phase, tuple()),
            candidate_evaluations=evaluations,
            crc_solution_count=1,
            repair_cost=selected.cost,
            frame=selected.frame,
            selected_phase=candidate.phase,
            selected_offset=candidate.offset,
            phases_tested=phases,
            beam_candidate_count=len(all_candidates),
            beam_flip_indices=candidate.flip_indices,
            syndrome_flip_destuffed_indices=selected.syndrome_destuffed_indices,
            syndrome_flip_raw_indices=selected.syndrome_raw_indices,
            syndrome_flip_count=len(selected.syndrome_destuffed_indices),
            arbitration_distinct_frames=distinct_count,
            arbitration_score=consensus_scores[winning_frame],
            arbitration_margin=arbitration_margin,
            repair_method="+".join(method_parts),
        )

    fallback = min(
        all_candidates,
        key=lambda candidate: (
            0 if not candidate.flip_indices else 1,
            candidate.cost,
            len(candidate.flip_indices),
            candidate.phase,
        ),
        default=None,
    )
    if competitive_count > 1:
        status = "AMBIGUOUS_CRC_SOLUTIONS"
        reason = (
            f"{competitive_count} minimum-weight CRC-valid frames survived "
            f"phase/path arbitration without the required soft-score margin "
            f"of {minimum_arbitration_margin:g} "
            f"({distinct_count} distinct CRC-valid frames total)."
        )
    elif no_flip_candidates:
        status = "STUFFING_ALIGNED_CRC_FAIL"
        reason = (
            "At least one phase has a strict 2192-bit HDLC path, but no unique "
            "1-4 bit syndrome repair passed CRC."
        )
    else:
        status = "UNRECOVERABLE"
        reason = (
            "No unique CRC-valid frame survived multi-phase HDLC beam search "
            "and 1-4 bit syndrome repair."
        )

    safe_fallback = fallback if fallback is not None and not fallback.flip_indices else None
    selected_phase = fallback.phase if fallback is not None else phase_items[0].phase
    selected_offset = fallback.offset if fallback is not None else phase_items[0].offset
    return AdvancedRepairResult(
        status=status,
        reason=reason,
        original_stuffing_aligned=original_aligned,
        original_crc_pass=original_crc,
        hard_anchor_match=bool(fallback),
        closing_flag_found=bool(fallback),
        selected_flag_start=fallback.flag_start if fallback is not None else None,
        destuffed_bit_count=EXPECTED_FRAME_BITS if fallback is not None else 0,
        excluded_indices=(
            safe_fallback.excluded_indices if safe_fallback is not None else tuple()
        ),
        flip_indices=tuple(),
        source_indices=(
            safe_fallback.source_indices if safe_fallback is not None else tuple()
        ),
        candidate_pool=raw_pools.get(selected_phase, tuple()),
        candidate_evaluations=evaluations,
        crc_solution_count=competitive_count,
        repair_cost=None,
        frame=safe_fallback.frame if safe_fallback is not None else None,
        selected_phase=selected_phase,
        selected_offset=selected_offset,
        phases_tested=phases,
        beam_candidate_count=len(all_candidates),
        beam_flip_indices=tuple(),
        syndrome_flip_destuffed_indices=tuple(),
        syndrome_flip_raw_indices=tuple(),
        syndrome_flip_count=0,
        arbitration_distinct_frames=distinct_count,
        arbitration_score=(
            consensus_scores[ranked_frames[0]] if ranked_frames else None
        ),
        arbitration_margin=arbitration_margin,
        repair_method="MULTIPHASE_HDLC_BEAM+CRC_SYNDROME_1TO4+UNIQUE_ARBITRATION",
    )
