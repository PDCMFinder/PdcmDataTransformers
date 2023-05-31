import requests
import json

class API:
    def __init__(self):
        link = ''

    def request(self):
        response_API = requests.get(self.link)
        return response_API
