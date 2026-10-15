# SCIONX CS16 MATLAB 可攜解碼工具包

## 解壓縮後如何執行

1. **先完整解壓 ZIP**，保留 `SCIONX_Decode_Toolkit` 內的目錄結構，不能只取出單支 `.m`。
2. 使用 MATLAB R2026a 開啟 `SCIONX/封包頻譜分析/Decode_Process/run_decode_process.m`。
3. **和原機相同：直接在此腳本修改 `inputDirectory`，儲存後按 Run。** 也可將 Current Folder 切到該 `Decode_Process` 目錄，執行 `run_decode_process`。預設處理隨包附上的 `20260720_220206_100000_98266_x.cs16`，並匯出圖表。
4. 結果寫在 `SCIONX/封包頻譜分析/Decode_Process/Results/Run_.../`；MATLAB workspace 的 `decodeResult` 是本次結果。

工具包可放在其他磁碟或目錄。執行目錄須可寫入，MATLAB 需正常啟用授權及 Java。原始輸入合計約 856 MiB，請預留至少 2 GB 空間供解壓與輸出。已驗證版本為 MATLAB R2026a；未逐一驗證其他 MATLAB 版本。

## 在另一台電腦檢查是否得到相同解碼

將 Current Folder 切到工具包根目錄，執行：

```matlab
report = selftest_portable_toolkit;
```

測試會先核對全部隨附檔案的大小與 SHA-256，再分別解碼 7/20 及 10/1 的錄音。每份錄音僅作單一比較輸入，7/23 Reference 由程式自動加入。測試逐一核對 burst 數、狀態、長度、CRC、符元率、停止原因、錨點時間（容差 1 μs），以及**還原出的每個 byte**是否與原機基準相同。成功時顯示 `PASS`，報告存為 `Verification/Portable_Test_Report.json`。

若只需先確認 ZIP 傳輸完整，可執行 `verify_toolkit_integrity`。這些測試會建立新的 Run 並更新 `Latest_Result.mat`。若修改了套件內的程式，完整性檢查會提示 SHA-256 不符；可以保留一份未修改的解壓副本作比對。

| 隨附錄音 | 用途／原機基準 |
| --- | --- |
| `20260723_091639_100000_98266_x.cs16` | **固定 Reference，完整附上，308658176 bytes**；參考控制為 274 bytes、CRC PASS。 |
| `20260720_220206_100000_98266_x.cs16` | 入口預設測試；1 burst、238 bytes、`PARTIAL_DIAGNOSTIC`、`UNVERIFIABLE_PARTIAL`。 |
| `20261001_221509_100000_69880.cs16` | 先前指定測試；1 burst、159 bytes、`PARTIAL_DIAGNOSTIC`、`UNVERIFIABLE_PARTIAL`。 |

部分封包的預期結果本來就是診斷輸出；還原內容與原機相同，不表示取得完整 frame 或 CRC PASS。原機比較基準保存在 `EXPECTED_RESULTS.json`，只包含比較錄音的候選 bytes，不匯出 Reference 的 TLM 或 bitstream。

## 改用自己的錄音或資料夾

編輯 `SCIONX/封包頻譜分析/Decode_Process/run_decode_process.m` 的 `inputDirectory`：

```matlab
inputDirectory="E:\captures";                 % 資料夾
inputDirectory="E:\captures\recording.cs16";  % 單一檔案
inputDirectory="";                            % 開啟資料夾選擇器
```

`recursive=true` 包含子資料夾，單檔模式會忽略它。`sampleRate=100000` 強制使用 100 kS/s；不同取樣率應設定實際數值，或設 `NaN` 讓程式從檔名推斷。`exportPlots=false` 可省略圖表輸出。Reference 的位置與取樣率由工具包內相對路徑固定設定，不受比較錄音的取樣率選項影響。

也可以直接呼叫：

```matlab
addpath(fullfile(pwd,'SCIONX','封包頻譜分析','Decode_Process'));
r = decode_scionx_folder("E:\captures\recording.cs16",SampleRate=100000);
```

## 流程與檔案

`run_decode_process.m` → `decode_scionx_folder` → `analyze_scionx_packets`。後者掃描 RF bursts、定位 postamble closing flags、分析重複資料、從錨點向前還原 HDLC、檢查 CRC，最後輸出 TLM。重複資料一致率只作診斷，不是這個入口的解碼門檻。正常解碼只需操作 `run_decode_process.m`；根目錄的完整性與比對測試是額外檢查工具。

完整且為 274-byte 的 frame 輸出 `TLM.txt`，CRC FAIL 也會保留；部分／較短 frame 輸出 `TLM_DIAGNOSTIC.txt`。先看 `Recording_Summary.csv` 及 `Packet_Summary.csv`，再依其中的路徑查看文字、對齊圖與診斷。

```text
SCIONX_Decode_Toolkit/
  README_zh-TW.md
  selftest_portable_toolkit.m
  verify_toolkit_integrity.m
  EXPECTED_RESULTS.json
  PACKAGE_MANIFEST.json
  SCIONX/
    Beacon Data/CS16_Raw/             # Reference 與兩份測試錄音，依日期分層
    SignalSim/GNU_Matlab_Demo/        # decode_scionx_tlm.m
    封包頻譜分析/                     # 共用 MATLAB 訊號處理函式
      Decode_Process/
        run_decode_process.m         # 直接修改 inputDirectory，然後執行
        decode_scionx_folder.m       # 資料夾／單檔包裝
        Results/                     # 執行後產生
```

CS16 解碼路徑使用 MATLAB 自帶功能，不需安裝 Node.js、Codex、Python 或額外 MATLAB toolbox。共享分析程式仍含其他輸入格式的分支；本工具包的範圍與測試對象為 CS16。

## 打包時的調整

核心演算法與 `decode_scionx_folder.m`、`decode_scionx_tlm.m` 直接複製原程式。工具包仍以 `run_decode_process.m` 為執行入口，僅將隨附錄音的預設位置改由自身路徑組合，並保留打包時的遞迴、取樣率、繪圖及重複 postamble 選項；另附完整性檢查及跨位置比對測試。原工作目錄的程式及既有 Results 不因打包而改寫。
