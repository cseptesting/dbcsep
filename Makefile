SOURCES=dbcsep tests setup.py
LENGTH=96

check: $(SOURCES)
	flake8 --max-line-length=$(LENGTH) $^
	black --check --line-length $(LENGTH) $^
	pylint -E $^

format: $(SOURCES)
	black --line-length $(LENGTH) $^
