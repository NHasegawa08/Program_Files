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
- (追記20220704) R版の関数identifyの挙動が良くない(図をクリックしても点が選択できない)ことなどから、python版を作りました。(pre-edge_analysis.ipynb)

## 2022/1/13 update XANES_spectrum_Plot.ipynb ~~& XANESFunc.py~~
- for python only
- BL-12Cのスペクトル描画用です。~~(2ファイルで１つ)~~
- データの個数に関わらず好きな数だけスペクトルを描けます。
- (20220708)内容を改変しました。
## 2022/6/8 XANES Fitting
- XANESスペクトルの~~2-3成分~~の線形結合フィッティングプログラム
- 最適な組み合わせ自動探索機能付
- (20220710)任意の数の成分に対応できるよう改良しました。
## 2022/6/9 BL-15A_Mapping
- 15Aのヒートマップ描画のメモ
