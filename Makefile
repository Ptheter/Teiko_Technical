.PHONY: setup pipeline dashboard

setup:
	pip install -r requirements.txt

pipeline:
	python3 clean_data.py
	python3 load_data.py
	python3 analysis.py

dashboard:
	streamlit run dashboard.py