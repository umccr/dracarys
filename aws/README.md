# dracarys CDK app

- [dracarys CDK app](#dracarys-cdk-app)
  - [Development](#development)
    - [Cfn Outputs](#cfn-outputs)
    - [Deployment Parameters](#deployment-parameters)
    - [Commands](#commands)
    - [Files](#files)
    - [Docker](#docker)

## Development

### Cfn Outputs

e.g. retrieve name of bucket

```typescript
const buck1 = new Bucket(this, 'aBucket');
new CfnOutput(this, 'myBucket', {
  value: buck1.bucketName,
});
```

Then after `cdk deploy` you have in the output log:

```text
Outputs:
MyCdkStack.myBucket = foo
```

### Deployment Parameters

e.g. to change the S3 bucket object expiration date on deploy:

```typescript
const duration_param = new CfnParameter(this, 'duration', {
  type: 'Number',
  default: 6,
  minValue: 1,
  maxValue: 10,
});

const buck2 = new Bucket(this, 'aBucket', {
  lifecycleRules: [
    {
      expiration: Duration.days(duration_param.valueAsNumber),
    },
  ],
});
```

### Commands

| Command                              | Description                                                        |
| ------------------------------------ | ------------------------------------------------------------------ |
| `cdk init app --language typescript` | initialise app from scratch (pretty much wraps the above commands) |
| `npm i -g aws-cdk`                   | installation of CDK globally                                       |
| `npm run build`                      | compile typescript to js                                           |
| `npm run watch`                      | watch for changes and compile                                      |
| `npm run test`                       | perform the jest unit tests                                        |
| `cdk deploy`                         | deploy this stack to your default AWS account/region               |
| `cdk diff`                           | compare deployed stack with current state                          |
| `cdk synth`                          | emits the synthesized CloudFormation template                      |

### Files

| File/Dir          | Description |
| ----------------- | ----------- |
| `cdk.json`        | execute app |
| `bin/dracarys.ts` | launch app  |
| `lib/`            | CDK stack   |

### Docker

Interactive debugging:

```bash
docker run --rm -it \
  --entrypoint /bin/bash \
  --platform linux/amd64 \
  public.ecr.aws/lambda/python:3.9
```