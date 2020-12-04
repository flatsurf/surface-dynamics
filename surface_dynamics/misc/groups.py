def relations(x, y):
    r"""
    Iterate through the relations between ``x`` and ``y`` by length
    """
    lx = 'x'
    ly = 'y'
    lX = 'X'
    lY = 'Y'
    if x.is_one():
        return lx
    if y.is_one():
        return ly
    gens = [(lx, x), (ly, y), (lX, ~x), (lY, ~y)]
    inv = {lx: lX, lX:lx, ly:lY, lY:ly}
    level = gens
    while True:
        new_level = []
        for w, g in level:
            for a, s in gens:
                if w[-1] != inv[a]:
                    gg = g * s
                    ww = w + a
                    if gg.is_one():
                        yield ww
                    new_level.append((ww, gg))
        level = new_level
