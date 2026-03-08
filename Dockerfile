FROM python:3.12-slim
WORKDIR /app
COPY . .
RUN pip install --no-cache-dir -r requirements.txt
# TODO: specify your entrypoint, e.g.:
# CMD ["python", "main.py"]
