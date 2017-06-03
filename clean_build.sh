git reset --hard HEAD

rm src/*.c *.so -rf build

python setup_cython.py build
python setupegg.py develop
python setup.py build
