.PHONY: all lint test test-cov install dev

all: 

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	py.test --cov=q2_ebd

q2_diversity/_alpha/alpha_correlation_assets/dist:
	cd q2_diversity/_alpha/alpha_correlation_assets && \
	npm install && \
	npm run build && \
	cp licenses/* dist/

q2_diversity/_alpha/alpha_group_significance_assets/dist:
	cd q2_diversity/_alpha/alpha_group_significance_assets && \
	npm install && \
	npm run build && \
	cp licenses/* dist/

q2_diversity/_alpha/alpha_rarefaction_assets/dist:
	cd q2_diversity/_alpha/alpha_rarefaction_assets && \
	npm install && \
	npm run build && \
	cp licenses/* dist/

install: all
	python setup.py install

dev: all
	pip install -e .
