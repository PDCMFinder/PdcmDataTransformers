from CancerModelsFinder.classes import API
from CancerModelsFinder.provider.PDXNet import transform
class PDXNet():
    def __init__(self):
        self.API = API.API()
        self.transform = transform.PDXNet_data_transformations()

    def request(self):
        response = self.API.request()

    def main(self):
        self.transform.process_models()
        self.transform.compare_models()
        if self.transform.flag:
            self.transform.transform()
