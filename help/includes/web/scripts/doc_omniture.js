// Copyright 2013 The MathWorks, Inc.
var TocTracking = {};
(function () {
    TocTracking.open = isTocOpen();
    TocTracking.clicked = isTocClicked();
    setCookie('MW_toc_clicked', false);

    function isTocOpen() {
        return getCookie('MW_toc_docked') === 'true';
    }

    function isTocClicked() {
        return getCookie('MW_toc_clicked') === 'true';
    }

    function getCookie(name) {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = cookies[i].replace(/^\s+|\s+$/g, '');
            if (cookie.indexOf(name) === 0) {
                return cookie.substring(name.length + 1, cookie.length);
            }

        }
        return null;
    }

    function setCookie(name, value) {
        var date = new Date();
        date.setTime(date.getTime() + (7 * 24 * 60 * 60 * 1000));
        var expiresDate = date.toGMTString();
        document.cookie = name + "=" + value
            + "; expires=" + expiresDate
            + "; path=/";

    }

    $(document).on('toc_clicked', function() {
        setCookie('MW_toc_clicked', true);
    });
})();

