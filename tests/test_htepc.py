import HTEPC as ht

obj = ht.MP_connect()
x = obj.setting("mp-763")
assert x == "mp-763"
