# SigmaP
Python package for Sigma70 promoter Prediction. This package used Sigma70Pred [(Patiyal et al. 2022)](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2022.1042127/full).

### Installation
This package can be installed by pip.
```python
pip install sigmap
```

### How to use
First, prepare fasta file containing DNA sequence. Minimum length for prediction is 81nt. Then, calculate probability score by `SigmaFactor`. Run prediction model by `.predict` method. Results will be returned as `pd.DataFrame`.
```python
from sigmap import SigmaFactor

sigma = SigmaFactor()

df_out = sigma.predict('tutorial/example_seq.fa')
```

| ID      | Sequence                                          | Score             | Prediction   |
| ------- | ------------------------------------------------- | ----------------- | ------------ |
| \>Seq_1 | TAGCACGACGATAATATAAACGCAGCAAAAAAAAAAAAAAAAAAAA... | 0.145             | Non-Promoter |
| \>Seq_2 | AGCTTGCGTCAATGGGCAAGGTGGGCTTGCATTTGCTTAATAGAAA... | 0.478             | Promoter     |
| \>Seq_3 | TCGTTTTATTTCTTTTTTCTCCATTGAACTTTCAGTTTCTTTTCTA... | 0.692             | Promoter     |
| \>Seq_4 | CGCAGCGGGTTTACCCTCTGACCGTTTCTGTTACGAAGGCTTTTTA... | 0.216             | Non-Promoter |
| \>Seq_5 | TGCTGCTTGGTCTGTGGGTTGCCGCACAGGTTGCCGGTTCCACCAA... | 0.162             | Non-Promoter |
| \>Seq_6 | GAATCCAACTAATGTTGTAAACTGGCAAGGTAATGTCATTAGTCAT... | 0.418             | Promoter     |


The input type for sigmap can also be a `pd.DataFrame`. If you want to convert a FASTA file into a DataFrame, you can use the `fasta2df` function.

```python
from sigmap import fasta2df

df_seq = fasta2df('tutorial/example_seq.fa')
```

| Sequence_ID | Sequence                                          |
| ----------- | ------------------------------------------------- |
| \>Seq_1     | TAGCACGACGATAATATAAACGCAGCAA                      |
| \>Seq_2     | AGCTTGCGTCAATGGGCAAGGTGGGCTTGCATTTGCTTAATAGAAA... |
| \>Seq_3     | TCGTTTTATTTCTTTTTTCTCCATTGAACTTTCAGTTTCTTTTCTA... |
| \>Seq_4     | CGCAGCGGGTTTACCCTCTGACCGTTTCTGTTACGAAGGCTTTTTA... |
| \>Seq_5     | TGCTGCTTGGTCTGTGGGTTGCCGCACAGGTTGCCGGTTCCACCAA... |
| \>Seq_6     | GAATCCAACTAATGTTGTAAACTGGCAAGGTAATGTCATTAGTCAT... |


If the `DataFrame` contains data with ID and sequence columns, you can directly use it as input for `SigmaFactor`.

```python
sigma = SigmaFactor()

# input type: pd.DataFrame
df_out = sigma.predict(df_seq)
```

| ID      | Sequence                                          | Score             | Prediction   |
| ------- | ------------------------------------------------- | ----------------- | ------------ |
| \>Seq_1 | TAGCACGACGATAATATAAACGCAGCAAAAAAAAAAAAAAAAAAAA... | 0.145             | Non-Promoter |
| \>Seq_2 | AGCTTGCGTCAATGGGCAAGGTGGGCTTGCATTTGCTTAATAGAAA... | 0.478             | Promoter     |
| \>Seq_3 | TCGTTTTATTTCTTTTTTCTCCATTGAACTTTCAGTTTCTTTTCTA... | 0.692             | Promoter     |
| \>Seq_4 | CGCAGCGGGTTTACCCTCTGACCGTTTCTGTTACGAAGGCTTTTTA... | 0.216             | Non-Promoter |
| \>Seq_5 | TGCTGCTTGGTCTGTGGGTTGCCGCACAGGTTGCCGGTTCCACCAA... | 0.162             | Non-Promoter |
| \>Seq_6 | GAATCCAACTAATGTTGTAAACTGGCAAGGTAATGTCATTAGTCAT... | 0.418             | Promoter     |



Contact: Goosang Yu (gsyu93@gmail.com)