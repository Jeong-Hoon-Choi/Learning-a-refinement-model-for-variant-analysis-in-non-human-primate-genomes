# Learning a Refinement Modle for Variant Analysis in Non-Human Primate Genomes

This repository contains the source code used to train and evaluate the genome variant refinement model proposed in the paper:  
**"Learning a Refinement Model for Variant Analysis in Non-Human Primate Genomes"**

It supports model development, hyperparameter tuning, and evaluation using both human and non-human primate genome datasets.

---

## 📁 Repository Structure

- `src/` — Source code of all the experiments that contains the training & evluation  
- `model/` — Training code for various machine learning models including XGBoost, LightGBM, and other classifiers.  

---

## 🔧 Dependencies

This project uses:

- Python 3.8+
- `xgboost`, `lightgbm`, `scikit-learn`, `pysam`, `seaborn`, `matplotlib`, `numpy`, `pandas`, etc.

Install all dependencies with:

```bash
pip install -r requirements.txt
```

---

## 🧪 Supported Experiments

This repository supports:

- Cross-sample evaluation (e.g., training on HG001, testing on HG002)
- Mixed training with HG001 and HG002
- Evaluation on Rhesus Macaque individuals
- Feature ablation studies
- Alternative base ratio (ABR) analysis and visualization
- Comparative benchmarking with other ML/DL models (e.g., RF, LR, k-NN, MLP, FT-Transformer)
