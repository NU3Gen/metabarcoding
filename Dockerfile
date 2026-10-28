FROM python:3.12-slim
WORKDIR /app
COPY pyproject.toml README.md ./
COPY src ./src
COPY references ./references
RUN pip install --no-cache-dir .
RUN useradd --create-home appuser && mkdir -p /app/results && chown -R appuser:appuser /app
USER appuser
EXPOSE 8000
CMD ["uvicorn", "metabarcoder.webapp:app", "--host", "0.0.0.0", "--port", "8000"]

