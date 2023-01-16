#!/usr/bin/env node
import 'source-map-support/register';
import * as cdk from 'aws-cdk-lib';
import { DracarysStack } from '../lib/dracarys-stack';

const app = new cdk.App();
new DracarysStack(app, 'DracarysStack', {
  description: 'Stack for dracarys',
  tags: {
    ['developer']: 'Peter Diakumis',
    ['project']: 'dracarys',
  },
});
