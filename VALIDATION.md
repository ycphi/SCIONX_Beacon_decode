# 驗證記錄

測試日期：2026-09-14。本機 Windows x64、MATLAB R2026a、內附 CPython 3.14.6。
MATLAB R2024a/R2024b 為預期目標；本機沒有這兩版，因此未宣稱已實機驗證。

- `scionx_portable_scan.py --check-environment` 通過：NumPy 2.4.4、SciPy 1.17.1、SoundFile 0.13.1、libsndfile 1.2.2。
- `gnuradio` 未安裝於內附環境，也未載入；soft repair 從套件內 `soft_repair` 載入。
- CRC-16/X-25 標準向量 `123456789` 結果為 `0x906E`。
- 封裝核心 Python、兩個 MATLAB decoder 與 JSON 的 SHA256 和本機版本完全一致。
- `run_scionx_pipeline` 僅調整內附 runtime 優先順序與相關說明／錯誤提示。
- 以 MATLAB 入口同時輸入 CS16、OGG，process pool 使用 2 workers，兩個 Python exit code 均為 0。

| 測試片段 | Frame | CRC PASS | TLM | 結果 |
|---|---:|---:|---:|---|
| 20260723_091639_100000_98266_x.cs16，82–87 秒 | 1 | 1 | 1 | SUCCESS，無需修復 |
| satnogs_14459039_2026-07-07T09-57-46.ogg，69–73 秒 | 1 | 1 | 1 | SUCCESS_WITH_SOFT_REPAIR |

CS16 保留原始 IQ；OGG 片段經 Vorbis 重新編碼，不能用本次結果推論完整原始 OGG 的修復位元數。
測試資料與輸出未納入 ZIP。這是流程驗證，不代表所有錄音皆能解碼或修復結果必然正確。

ZIP 內 `MANIFEST_SHA256.json` 保存所有封裝檔案（不含 manifest 本身）的雜湊；建置時逐檔核對 ZIP 內容。
