from multiCompApp.ua import ua_client


def get_all_rxs():
    ua = ua_client.UaClient()
    while True:
        data = ua.get_data(0)  # zależy w jakim slocie to było zapisane
        reactions = [x.decode('utf-8') for x in data.split(b'\n')]
        if len(reactions) > 2:
            break
        time.sleep(2)
    ua.close()
    return reactions


rxadb = get_all_rxs()