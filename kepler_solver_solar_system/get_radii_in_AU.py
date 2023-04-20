from astropy import units as u
#venus, mars, saturn, uranus, neptune
d_in_km = [12104,6792,120536,51118,49528]

for d in d_in_km:
    r = d/2.0
    r_in_km = r * u.km
    r_in_au = r_in_km.to(u.AU).value
    print("{:21.16f}".format(r_in_au))

