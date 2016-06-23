echo 'Running py.test...'
py.test --cov --cov-report term-missing

echo ''
echo 'Running pycodestyle test...'
pycodestyle .

echo '----'
echo 'End test.'