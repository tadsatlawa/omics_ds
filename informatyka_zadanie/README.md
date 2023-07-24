# Rozwiązanie zadania

1. `./get_files.sh`
2. `./extract_sequences.sh mouse.1.protein.faa sample_100 ">" 100`
3. `./create_db.sh zebrafish.1.protein.faa zebrafish_db`
4. `sbatch -p topola script.sbatch` --> creates out file and slurm-*.out
5. `./count_sequences.sh mouse.*`  --> 90730

Output dla próbki 100 sekwencji
- rozmiar: 32MB
-
```
real	0m17.970s  # 24 threads
real	0m53.998s  # 4 threads
real	2m59.828s  # 1 treads
``````

# Odpowiedzi na pytania

1. Obliczyć czas (w godzinach) potrzebny do	wykonania obliczeń dla pełnej sekwencji RNA
myszy.

- przy załozeniu uruchomienia na 24 wątkach czas wykonania, to ~4.5h

2. Obliczyć	przestrzeń dyskową	potrzebną do zapisu	wyników	(w GB).

- ~29 GB

3. Ile razy	szybciej (w	stosunku do wykorzystania 1 rdzenia obliczeniowego) wykonamy obliczenia korzystając z 4 rdzeni obliczeniowych, a ile razy szybciej korzystając z 24 rdzeni obliczeniowych.

- 1 -> 4: 3.4x
- 1 -> 24: 10x
