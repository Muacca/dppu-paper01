# DPPU Paper 1: Microscopic Handles in Teleparallel Gravity（日本語）

このファイルは日本語版の README です。  
英語版 README は以下をご覧ください：  
➡️ [English README](./README.md)

---

## 概要

このリポジトリは、DPPU プロジェクトにおける第一論文  
「Microscopic Handles in Teleparallel Gravity」の

- 論文本体の LaTeX ソース（本文＋付録）
- Appendix E で用いた SymPy スクリプトとそのログ

を公開するためのものです。

第一論文では、テレパラレル重力（TEGR）の枠組みの中で、トーショナル単極電荷をもつ微視的ハンドル（$S^2 \times S^1$ ）を用いて、

- 安定なハンドル半径と $E \propto |q|$ な質量スケール
- Nieh–Yan 項を介したカイラリティ（左右非対称性）の選択
- 前進モードの剛性 $k(q) \propto \omega^2 |q|$ と古典的臨界指数 $\gamma = 1$

といった性質が幾何学的にどこまで説明できるかを検討しています。

---

## 構成

- `paper/`  
  論文の LaTeX ソース一式です（本文セクションと各 Appendix、参考文献、図ファイルなど）。

- `code/`  
  Appendix E で用いた SymPy スクリプトと、その実行ログです。  
  フェーズ 3 における剛性 $k(q)$ の導出や、$q^2 \varepsilon^2$ 項のキャンセル確認に対応しています。

---

## ビルド方法

LaTeX のビルドには `latexmk` の利用を想定しています。  
典型的なビルド手順は次の通りです：

```bash
cd paper
latexmk -pdf main.tex
````

依存パッケージ（`amsmath`, `graphicx`, `hyperref` など）が
TeX 環境に入っていることを前提としています。

---

## コードとデータの公開について

Appendix E の計算に対応する SymPy スクリプトと、
その実行結果のログはすべて `code/` ディレクトリに含まれています。

これには、たとえば以下のような処理が含まれます：

* トーションテンソルの摂動展開
* $q^2 \varepsilon^2$ 項の係数の自動抽出とキャンセル確認
* 剛性 $k(q)$ の $\omega^2 |q|$ スケーリングの検証 など

リポジトリの正規 URL は次の通りです：

`https://github.com/Muacca/dppu-paper01`

---

## ライセンス

このリポジトリでは、コードとテキストでライセンスを分けています：

* `code/` 配下の **すべてのコード** は MIT License の下で公開されています
  （詳細は `LICENSE-code` を参照してください）。
* `paper/` 配下の **テキスト・図・LaTeX ソース** は
  Creative Commons Attribution 4.0 International（CC BY 4.0）ライセンスの下で公開されています
  （詳細は `LICENSE-text` を参照してください）。


