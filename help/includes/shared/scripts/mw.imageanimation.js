(function($) {
    $.fn.animateImage = function() {

        var checkAnimatedImage = function (imgElt) {
            imgElt.data('src', imgElt.attr('src'));
            imgElt.click(toggleImageAnimation);
        };

        function toggleImageAnimation() {
            var imgElt = $(this);
            if (imgElt.hasClass('dynamic-image')) {
                //swap in the static image
                imgElt.attr('src', imgElt.data('src'));
                imgElt.removeClass('dynamic-image');
                imgElt.attr('title', getLocalizedString('play'));
            } else {
                //swap in the animated image.
                var anm_prefix = "_anmtd_";
                var anm_suffix = ".gif";
                var pathArray = imgElt.data('src').split("/");
                var fileName = pathArray.slice(-1)[0];
                pathArray[pathArray.length - 1] = anm_prefix + fileName.split(".")[0] + anm_suffix;
                var anm_src = pathArray.join("/");
                imgElt.attr('src', anm_src);
                imgElt.addClass('dynamic-image');
                imgElt.attr('title', getLocalizedString('stop'));
            }
        }

        return this.each(function () {
            checkAnimatedImage($(this));
        });

    };
})(jQuery);