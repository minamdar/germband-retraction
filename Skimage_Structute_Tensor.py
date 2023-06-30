# skimage structure tensor

def dominant_direction(img, sigma):
    """OrientationsJ's dominant direction"""
    axx, axy, ayy = feature.structure_tensor(
        img.astype(numpy.float32), sigma=sigma, mode="reflect"
    )
    dom_ori = numpy.arctan2(2 * axy.mean(), (ayy.mean() - axx.mean())) / 2
    return numpy.rad2deg(dom_ori)


def orientation_analysis(img, sigma):
    """OrientationJ's output for
    * orientation
    * coherence
    * energy
    """
    eps = 1e-20

    axx, axy, ayy = feature.structure_tensor(
        img.astype(numpy.float32), sigma=sigma, mode="reflect"
    )
    l1, l2 = feature.structure_tensor_eigvals(axx, axy, ayy)
    ori = numpy.arctan2(2 * axy, (ayy - axx)) / 2

    coh = ((l2 - l1) / (l2 + l1 + eps)) ** 2
    ene = numpy.sqrt(axx + ayy)
    ene /= ene.max()

    return ori, coh, ene