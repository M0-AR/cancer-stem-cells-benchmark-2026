run:
	python experiments/run_benchmark.py

compose:
	docker compose up --build

test:
	python -m py_compile src/csc_benchmark/*.py experiments/run_benchmark.py
