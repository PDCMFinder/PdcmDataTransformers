import requests

def download_file(url, file_name):
    """
    Download a file from a URL and save it to a local file.
    """
    with open(file_name, 'wb') as f:
        response = requests.get(url, stream=True)
        total_length = response.headers.get('content-length')

        if total_length is None: # no content length header
            f.write(response.content)
        else:
            dl = 0
            total_length = int(total_length)
            for data in response.iter_content(chunk_size=4096):
                dl += len(data)
                f.write(data)
                done = int(50 * dl / total_length)
                print("\r[%s%s]" % ('=' * done, ' ' * (50-done)), end='')
            print()