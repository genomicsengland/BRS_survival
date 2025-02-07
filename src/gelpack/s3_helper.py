import boto3
import pandas as pd

def get_s3():
    annot_genome_access_role = (
        "arn:aws:iam::907999473992:role/Gel_Data_Resources_MM_S3_Read_Access"
    )
    session_name = "gel-data-access"

    sts = boto3.client("sts")
    response = sts.assume_role(
        RoleArn=annot_genome_access_role, RoleSessionName=session_name
    )
    new_session = boto3.Session(
        aws_access_key_id=response["Credentials"]["AccessKeyId"],
        aws_secret_access_key=response["Credentials"]["SecretAccessKey"],
        aws_session_token=response["Credentials"]["SessionToken"],
    )
    s3 = new_session.client("s3")

    return s3
    
def load_vcf_s3_to_local(file_path, local_path):
    """
    Download a VCF file from S3 to a local path.

    Parameters:
        file_path (str): The path of the file in the S3 bucket.
        local_path (str): The local path where the file will be downloaded.
    """
    
    s3 = get_s3()

    bucket = file_path.replace("s3://", "").split("/")[0]
    key = file_path.replace("s3://", "").replace(bucket, "")[1:]

    try:
        s3.download_file(Bucket=bucket, Key=key, Filename=local_path)
    except Exception as e:
        print(f"Error occurred while downloading the file: {e}")


def read_tsv_from_rdl(s3path, sample=False):
    role_to_assume_arn = "arn:aws:iam::907999473992:role/RDL_MM_S3_Read_Access"

    sts_client = boto3.client("sts")

    # Assume the role
    response = sts_client.assume_role(
        RoleArn=role_to_assume_arn, RoleSessionName="AssumeRoleSession"
    )

    # Extract the temporary credentials from the response
    credentials = response["Credentials"]

    # To create a boto3 session for download_file etc.
    s3_session_credentials_dict = {
        "key": credentials["AccessKeyId"],
        "secret": credentials["SecretAccessKey"],
        "token": credentials["SessionToken"],
    }

    if "s3://" not in s3path:
        bucket = '676167504555-gel-structured-data-prod/'
        s3path = 's3://'+bucket+s3path

    try:
        if sample:  # Check if sampling is enabled
            df = pd.read_csv(s3path, storage_options=s3_session_credentials_dict, sep="\t", nrows=100)
        else:  # Read the entire file
            df = pd.read_csv(s3path, storage_options=s3_session_credentials_dict, sep="\t")
    except UnicodeDecodeError:
        encoding = "latin-1"
        if sample:  # Apply the same logic for the exception case
            df = pd.read_csv(s3path, storage_options=s3_session_credentials_dict, sep="\t", encoding=encoding, nrows=100)
        else:
            df = pd.read_csv(s3path, storage_options=s3_session_credentials_dict, sep="\t", encoding=encoding)

    return df

def get_rdl_catalogue():
    role_to_assume_arn = 'arn:aws:iam::907999473992:role/RDL_MM_S3_Read_Access'
    
    # Create an STS client
    sts_client = boto3.client('sts')
    
    # Assume the role
    response = sts_client.assume_role(
    RoleArn=role_to_assume_arn,
    RoleSessionName='AssumeRoleSession'
    )
    
    # Extract the temporary credentials from the response
    credentials = response['Credentials']
    
    # To create a boto3 session for download_file etc.
    session_credentials_dict = {'aws_access_key_id': credentials['AccessKeyId'],
    'aws_secret_access_key': credentials['SecretAccessKey'],
    'aws_session_token': credentials['SessionToken']}
    
    # Create a new session using the assumed role credentials
    assumed_session = boto3.Session(**session_credentials_dict)
    
    # Now you can use the assumed session for any AWS service
    s3_assumed_client = assumed_session.client('s3')

    paginator = s3_assumed_client.get_paginator("list_objects_v2")
    pages = paginator.paginate(Bucket='676167504555-gel-structured-data-prod')
    conts = []
    for page in pages:
        for obj in page["Contents"]:
            conts.append(obj)
    rdl_catalogue = pd.DataFrame(conts)

    return rdl_catalogue