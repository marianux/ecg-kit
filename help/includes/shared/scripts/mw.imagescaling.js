(function($) {
    $.fn.scaleImage = function(options) {
        options = $.extend({}, defaults, options);

        var checkImage = function (imgElt) {
            var callback = function () {
                if (needsScaling(imgElt)) {
                    imgElt.addClass('scaled-image');
                    imgElt.click(toggleImageScaling);
                }
            };
            //Preload image.
            loadImage(imgElt, callback);
        };

        function needsScaling(imgElt) {
            var fullSize = imgElt.data('fullsize');
            //Vary width and height values here to change scaling parameters.
            return fullSize.width > options.width; /*|| fullSize.height > options.height;*/
        }

        function loadImage(imgElt, callback) {
            var img = new Image();
            img.onload = function () {
                imgElt.data('fullsize', {
                    'width':img.width,
                    'height':img.height
                });
                callback();
            };
            img.src = imgElt.attr('src');
        }

        //Setup image scaling.
        function toggleImageScaling() {
            var imgElt = $(this);
            if (imgElt.hasClass('expanded-image')) {
                imgElt.animate({
                    'max-width':options.scaledWidth,
                    'max-height': options.scaledHeight
                }, function () {
                    imgElt.removeClass('expanded-image');
                });
            } else {
                var fullSize = imgElt.data('fullsize');
                imgElt.animate({
                    'max-width':fullSize.width + 'px',
                    'max-height':fullSize.height + 'px'
                }, function () {
                    imgElt.addClass('expanded-image');
                });
            }
        }

        return this.each(function () {
            checkImage($(this));
        });
    };

    var defaults = {
        width: 600,
        height: 600,
        scaledWidth: '400px',
        scaledHeight: '300px'
    };
})(jQuery);