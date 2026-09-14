# SCIONX_Beacon_decode：Windows 可攜式信標解碼

## 開始使用

1. 將 ZIP **完整解壓縮**到可寫入的資料夾，例如 `D:\SCIONX_Beacon_decode`。不要只複製 `.m`，也不要在 ZIP 預覽中執行。
2. 在 MATLAB 2024a／2024b 的 Current Folder 切換到解壓縮後、含 `run_scionx_pipeline.m` 的資料夾。
3. 在 Command Window 執行：

```matlab
results = run_scionx_pipeline();
```

依對話框選擇一個 CS16／OGG 檔，或包含錄音檔的資料夾。程式會遞迴搜尋資料夾，可混合兩種格式。

也可以指定路徑：

```matlab
cd('D:\SCIONX_Beacon_decode');  % 改成實際解壓縮位置
results = run_scionx_pipeline('E:\Recordings');
results = run_scionx_pipeline('E:\Recordings\example.cs16', 'E:\Decode_Output');
% 指定 2 個平行 worker；第三個參數留空以使用內附 Python。
results = run_scionx_pipeline('E:\Recordings', 'E:\Decode_Output', [], 2);
% 使用序列處理：
results = run_scionx_pipeline('E:\Recordings', 'E:\Decode_Output', [], 1);
```

不需要原機的 D 槽路徑。預設輸出位於本套件的 `Pipeline_Output`。請為不同批次使用不同輸出資料夾，以免覆寫同名結果。

## 執行環境與相容性

- Windows 10／11 **x64（Intel／AMD）**，不是 32-bit 或 ARM64 原生套件。
- 需要已安裝並授權的 MATLAB；目標為 R2024a／R2024b。
- 內附獨立的 CPython **3.14.6 x64**、NumPy 2.4.4、SciPy 1.17.1、SoundFile 0.13.1、CFFI 2.0.0、pycparser 3.0 及音訊 DLL。
- 不需要 RadioConda、GNU Radio、pip 安裝或連網。另一台電腦原有 Python 3.14.6 不會被修改，也不必設定 PATH。
- 電腦原本安裝任何 Python 3 版本，甚至沒有 Python，都不影響預設執行；因為程式使用內附版本。這不代表程式碼可以交給任意舊版 Python 3 執行。
- pipeline 透過 `system` 啟動獨立 Python 程序，**不是 MATLAB 的 `pyenv`／`py.*` 介面**；不必將 Python 3.14 設為 MATLAB 的內嵌 Python。
- 平行處理需要 Parallel Computing Toolbox；沒有授權或無法建立 process pool 時，會退回序列處理。自動模式對資料夾最多使用 2 個 worker，MAXIMUM 搜尋可能耗時且耗記憶體。
- 此處實測的 MATLAB 版本請見 `VALIDATION.md`；未在您的另一台電腦實測，不能保證其作業系統／安全軟體配置。

## 程式分工與參數

`run_scionx_pipeline.m` → `scionx_portable_scan.py` → 共用 burst／候選／soft repair 核心 → `decode_gnuradio_hex.m` → `decode_scionx_tlm.m`。

Python 負責訊號解調、候選評分與 soft repair 證據；MATLAB 負責正式／診斷 HDLC deframe 與 TLM。診斷或修復後資料不代表無誤，請同時檢查 CRC 與修復狀態。

所有掃描參數修改 `decoder_settings.json`，修改後重新執行。沿用打包時本機設定：MAXIMUM soft repair、normal/inverted、9550–9650 sym/s、clock recovery、RAW/FIR。CS16 搜頻 ±20 kHz、100000 samples/s；OGG 搜頻 ±1000 Hz、48000 samples/s。實際值以 JSON 為準，不是所有錄音都適用。

- CS16：little-endian、交錯 signed int16 I/Q；sample rate 必須與 JSON 一致。
- OGG：單聲道解調音訊、48000 Hz；不是雙聲道 IQ。程式會拒絕與設定不符的聲道／取樣率。
- 不附原始錄音、歷史輸出、Excel 整合或 GRC 編輯器。這是可直接執行的 PORTABLE_PYTHON 版本，不提供 GNU_RADIO backend。
- 未改動訊號解調或遙測欄位規則；移植版僅將 runtime 選擇改為優先內附環境。NumPy/SciPy 與本機舊環境版本不同，邊界數值可能有差异，不承諾逐位元完全一致。

## 輸出與排錯

`Pipeline_Summary.csv` 記錄每個錄音的狀態、frame／CRC／TLM 數與沒有 TLM 的原因。每段錄音子資料夾下有 burst report、Matlab_Frame_Report.csv 與 `Frames/Frame_NNN`，可包含 HEX、Raw_bit.txt、TLM.txt 及 CRC／修復報告。沒有足夠候選時不一定產生 TLM。

如果 MATLAB 找錯版本，執行 `which run_scionx_pipeline -all`，確認指向解壓縮的檔案。第一次使用建議 `clear functions` 後重新呼叫。

可在套件資料夾的 PowerShell 檢查 Python：

```powershell
.\runtime\python.exe .\scionx_portable_scan.py --check-environment
```

預期出現 `soft_repair: AVAILABLE`、`gnuradio_imported: false`。若 runtime 遺失或 DLL 載入失敗，請重新完整解壓、確認防毒未隔離檔案、確認是 Windows x64。不要刪除 `runtime`、`soft_repair` 或更動其相對位置。

第三個參數可以明確指定其他 python.exe，但該環境必須自行安裝相容套件；免安裝保證只針對內附環境。若報錯，不要依舊版通用訊息安裝 RadioConda，先執行上述環境檢查。

## 第三方來源與授權

CPython 下載及校驗來源：https://www.python.org/downloads/release/python-3146/

NumPy／SciPy Python 3.14 支援說明：
https://numpy.org/doc/2.4/release/2.4.0-notes.html
https://docs.scipy.org/doc/scipy/release/1.17.0-notes.html

Python 授權保留於 `runtime/LICENSE.txt`；套件及所附 DLL 授權保留在 `runtime/Lib/site-packages` 的各套件與 `.dist-info` 目錄；soft repair 原始授權位於 `soft_repair/LICENSE`。
