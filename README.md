# Credit Risk & Credit Derivatives

**CDS curves · Bespoke indices · Tranche analysis · Bond ETF volatility**

Coursework from MF772 exploring how default risk enters credit-derivative pricing and how credit instruments behave in market data. The repository contains R and Python analyses with accompanying assignment reports.

## Explore the analyses

| Study | Topics | Code |
| --- | --- | --- |
| Bespoke CDX index | CDS survival-curve bootstrapping, premium and protection legs, tranche spreads, sensitivities, and jump-to-default calculations | [R analysis](Creat%20A%20Bespoke%20CDX%20Index/MF772%20HW3.R) |
| Credit indices & bond ETFs | Stylized hazard-rate calculations; HYG, JNK, and LQD return and rolling-volatility analysis | [Python notebook](Implied%20Volatility%20And%20the%20bespoke%20index/MF772%20HW2.ipynb) |

Read the [CDX assignment report](Creat%20A%20Bespoke%20CDX%20Index/mf772-hw3.pdf) or the [credit index and ETF report](Implied%20Volatility%20And%20the%20bespoke%20index/mf772-hw2.pdf) for the original questions and discussion.

## Methods & tools

- **R:** `credule` for CDS curve bootstrapping and `Ryacas0` for symbolic calculations.
- **Python:** pandas, NumPy, SymPy, and Matplotlib for data analysis and visualization.
- The R exercise examines how recovery and correlation assumptions relate to stylized portfolio and tranche calculations.

## Working with the files

The original coursework is preserved. The Python notebook contains local file paths that must be changed before execution. HYG, JNK, and LQD CSV files are included; a later section also references `SPY.csv`, which is not included in this repository. R package compatibility and input dimensions should be checked in the intended environment.

These are educational implementations with simplified assumptions, rather than independently validated pricing or risk systems. Results have not been rerun as part of this documentation update.

**Author:** Sike (Estelle) Yang
