## XANES_file_exchange
- This program combines multiple data files into one by extracting the necessary information and performing energy conversion from XANES raw file.
- Accept to BL-15A
## XANES_Data_Analyzer
- Accept to the data files of **XANES_file_exchange**
- This program conducts (1). removing abnormal points, (2). correcting BG, and (3). normalization of multiple raw data (15 defaults) at once and draw the normalized results. 
## XANES_pre-edge_Analyzer
- for pre-edge analysis
- fitting Gaussian function to the pre-edge of the XANES and determine the peak center.
- ***.xan*** file is available so far. 
## 2022/1/13 update XANES_spectrum_Plot.ipynb & XANESFunc.py
- for python only
- BL-12Cのスペクトル描画用です。(2ファイルで１つ)
- データの個数に関わらず好きな数だけスペクトルを描けます。
## 2022/6/8 XANES_Fitting
- XANESスペクトルの2-3成分の線形結合フィッティングプログラム
- 最適な成分の自動探索機能付
## 2022/6/9 BL-15A_Mapping
- 15Aのヒートマップ描画のメモ
