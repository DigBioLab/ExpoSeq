import onedrivesdk
from onedrivesdk.helpers import GetAuthCodeServer



redirect_uri = 'http://localhost:8080/'
client_secret = 'your_client_secret'
client_id='your_client_id'
api_base_url='https://api.onedrive.com/v1.0/'


scopes=['wl.signin', 'wl.offline_access', 'onedrive.readwrite']

client = onedrivesdk.get_default_client(
    client_id='your_client_id', scopes=scopes)

auth_url = client.auth_provider.get_auth_url(redirect_uri)

