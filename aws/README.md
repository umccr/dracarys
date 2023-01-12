# dracarys CDK app

## Development

### Commands

| Command                              | Description                                          |
| ------------------------------------ | ---------------------------------------------------- |
| `cdk init app --language typescript` | initialise app from scratch                          |
| `npm i -g aws-cdk`                   | installation of CDK globally                         |
| `npm run build`                      | compile typescript to js                             |
| `npm run watch`                      | watch for changes and compile                        |
| `npm run test`                       | perform the jest unit tests                          |
| `cdk deploy`                         | deploy this stack to your default AWS account/region |
| `cdk diff`                           | compare deployed stack with current state            |
| `cdk synth`                          | emits the synthesized CloudFormation template        |

### Files

| File/Dir   | Description        |
| ---------- | ------------------ |
| `cdk.json` | execute app        |
| `bin/`     | initialise app     |
| `lib/`     | contains CDK stack |
