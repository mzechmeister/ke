## CORDIC-like algorithm for KE in python

The python module [ke.py](ke.py) is vectorised and uses numpy-ctypes to wrap a C implementation.

Example:

```python
>>> import ke
>>> M = [1.0907, 2, 3.1, 4]
>>> ke._E(M, 1.)
array([1.99999818, 2.55419595, 3.12079558, 3.57764002])
```

The python module [ke_cordic.py](ke_cordic.py) hosts the demonstration codes from the paper and is standalone.


### Installation python
```bash
svn export https://github.com/mzechmeister/ke/trunk/python ke_cordic_like
cd ke_cordic_like
python setup.py develop --user
```

### Deinstallation

```bash
python setup.py develop -u --user
```
