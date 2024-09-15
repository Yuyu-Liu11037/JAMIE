# Enviroment
```
# JAMIE
git clone https://github.com/Yuyu-Liu11037/JAMIE.git
cd JAMIE
git config --global user.email "eliu11037@gmail.com"
git config --global user.name "Yuyu-Liu11037"
conda create -n JAMIE python=3.9 -y
source activate JAMIE
pip install gdown leidenalg
pip install -r requirements.txt
pip install -e .
```

# Citeseq
```
python citeseq/data.py
python citeseq/citeseq.py
python citeseq/eval.py
```

# Remark
Fast training and imputing, no need to save data.