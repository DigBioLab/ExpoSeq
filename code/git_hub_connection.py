import requests
from pprint import pprint

# github username
username = "x4nth055"
# url to request
url = f"https://api.github.com/users/{username}"
# make the request and return the json
user_data = requests.get(url).json()
# pretty print JSON data
pprint(user_data)


repo.create_file("test.txt", "commit message", "content of the file")
# delete that created file
contents = repo.get_contents("test.txt")
repo.delete_file(contents.path, "remove test.txt", contents.sha)