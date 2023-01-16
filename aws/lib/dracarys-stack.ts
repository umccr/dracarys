import * as path from 'path';
import * as cdk from 'aws-cdk-lib';
import { Construct } from 'constructs';
import { Duration } from 'aws-cdk-lib';
import * as s3 from 'aws-cdk-lib/aws-s3';
import * as lambda from 'aws-cdk-lib/aws-lambda';
import * as iam from 'aws-cdk-lib/aws-iam';
// import * as sqs from 'aws-cdk-lib/aws-sqs';

export class DracarysStack extends cdk.Stack {
  constructor(scope: Construct, id: string, props?: cdk.StackProps) {
    super(scope, id, props);

    new s3.Bucket(this, 'dracarysBucket', {
      bucketName: 'umccr-pd-dracarys-testing',
    });

    // add role and permissions for lambda
    const lambda_role = new iam.Role(this, 'My Role', {
      assumedBy: new iam.ServicePrincipal('lambda.amazonaws.com'),
    });

    lambda_role.addManagedPolicy(
      iam.ManagedPolicy.fromAwsManagedPolicyName(
        'service-role/AWSLambdaBasicExecutionRole'
      )
    );

    const lambda_duration_mins = 3;
    new lambda.DockerImageFunction(this, 'runDracarysFunction', {
      functionName: 'runDracarysFunction',
      description: 'dracarys lambda',
      code: lambda.DockerImageCode.fromImageAsset(
        path.join(__dirname, '../..')
      ),
      timeout: Duration.minutes(lambda_duration_mins),
      role: lambda_role,
    });
  }
}
