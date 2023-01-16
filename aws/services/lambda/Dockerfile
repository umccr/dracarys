FROM pdiakumis/dracarys:0.6.0
LABEL maintainer="https://github.com/pdiakumis"

# env vars from public.ecr.aws/lambda/python:3.9
ENV LAMBDA_RUNTIME_DIR=/var/runtime
ENV LAMBDA_TASK_ROOT=/var/task

# SHELL ["conda", "run", "-n", "dracarys_env", "/bin/bash", "-c"]

COPY requirements.txt ./
RUN pip install -r requirements.txt
COPY lambda_code.py ./
CMD ["lambda_code.handler"]
#CMD ["dracarys.R"]
