

def grey_color_func(
        word, font_size, position,
        orientation, random_state=None,
        **kwargs):
    """wordcloud colouring"""
    return "hsl(360, 0%, 0%)"


def clean_names(messyname):
    """clean out initials in authorship names"""
    cleanname = ' '.join([w for w in messyname.split() if len(w) > 1])
    return cleanname
