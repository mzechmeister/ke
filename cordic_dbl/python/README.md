## CORDIC-dbl algorithm for KE in python

The python module [ke_dbl_i.py](ke_dbl_i.py) is vectorised and uses numpy-ctypes to wrap a C implementation.

Example:

```python
>>> import ke_dbl_i
>>> M = [1.0907, 2, 3.1, 4]
>>> ke_dbl_i.i_E(M, 1.)
array([1.99999818, 2.55419595, 3.12079558, 3.57764002])
```

The python module [ke_cordic_dbl.py](ke_cordic_dbl.py) hosts the demonstration codes from the paper and is standalone.


### Installation python
```bash
svn export https://github.com/mzechmeister/ke/trunk/cordic_dbl/python ke_cordic_dbl
cd ke_cordic_dbl
python setup.py develop --user
```

### Deinstallation

```bash
python setup.py develop -u --user
```
