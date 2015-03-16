(function($) {
    $.fn.mwBreadcrumb = function() {
        var _private = {
            checkBreadcrumb: function (breadcrumb) {
                this.processBreadcrumb(breadcrumb);
                this.checkAndCollapseBreadCrumb(breadcrumb);
                var availableBreadcrumbWidth = breadcrumb.width();
                breadcrumb.data('width', availableBreadcrumbWidth);
            },

            processBreadcrumb: function (breadcrumb) {
                //Add classes to first and last elements in the breadcrumb.
                breadcrumb.find('li:first').addClass('breadcrumb_first');
                breadcrumb.find('li:eq(1)').addClass('breadcrumb_product');
                breadcrumb.find('li:last').addClass('breadcrumb_last');

                var children = breadcrumb.children();
                breadcrumb.append($('<div class="breadcrumbs_outer"></div>'));
                $(".breadcrumbs_outer").append($('<div class="breadcrumbs_inner"></div>'));
                $(".breadcrumbs_inner").append(children.remove());
            },

            collapseBreadcrumbElement: function(elem) {
                if (elem.hasClass('breadcrumbs_collapsed')) {
                    return;
                }
                //Save original width and text of breadcrumb
                var actualWidth = elem.width() + 7; //account for negative padding;
                elem.data({width: actualWidth});
                var innerSpan = elem.find('span');
                innerSpan.data({text: innerSpan.text()});

                //collapse breadcrumb
                elem.css({'width': '20px'});
                innerSpan.text("...");

                elem.addClass('breadcrumbs_collapsed');

                //add hover event handler for collapsed breadcrumb
                this.addBreadcrumbEventHandler(elem);
            },

            addBreadcrumbEventHandler: function(elem) {
                elem.bind('mouseenter.breadcrumb',
                    function() {
                        var span = $(this).find('span');
                        span.text(span.data('text'));
                        $(this).stop();
                        $(this).animate({
                            width: $(this).data('width')
                        });
                        return false;
                    }).bind('mouseleave.breadcrumb', function() {
                        var span = $(this).find('span');
                        $(this).stop();
                        $(this).animate({
                            width: 20
                        }, function () {
                            span.text("...");
                        });
                        return false;
                    });
            },

            checkAndCollapseBreadCrumb: function(breadcrumb) {
                var availableBreadcrumbWidth = breadcrumb.width();
                var collapsibleBreadcrumbChildren = $(".breadcrumbs_collapsed:last").nextAll();
                if (collapsibleBreadcrumbChildren.length === 0) {
                    collapsibleBreadcrumbChildren = $(".breadcrumbs_inner ul li:gt(1)");
                }

                var that = this;
                collapsibleBreadcrumbChildren.each(function (i, elem) {
                    var actualBreadcrumbWidth = 0;
                    breadcrumb.find('li').each(function() {
                        actualBreadcrumbWidth += $(this).width();
                    });
                    if (actualBreadcrumbWidth <= availableBreadcrumbWidth) {
                        return false;
                    }
                    that.collapseBreadcrumbElement($(elem));
                });
            },

            checkAndExpandBreadCrumb: function(breadcrumb) {
                var availableBreadcrumbWidth = breadcrumb.width();

                var expandableBreadcrumbChildren = $($(".breadcrumbs_collapsed").get().reverse());
                var that = this;
                expandableBreadcrumbChildren.each(function (i, elem) {
                    that.expandBreadcrumbElement($(elem));
                    var actualBreadcrumbWidth = 0;
                    breadcrumb.find('li').each(function() {
                        actualBreadcrumbWidth += $(this).width();
                    });
                    if (actualBreadcrumbWidth > availableBreadcrumbWidth) {
                        that.collapseBreadcrumbElement($(elem));
                        return false;
                    }
                });
            },

            expandBreadcrumbElement: function (elem) {
                //Save original width and text of breadcrumb
                var innerSpan = elem.find('span');
                innerSpan.text(innerSpan.data('text'));
                elem.css('width', 'auto');
                elem.removeClass('breadcrumbs_collapsed');
                elem.unbind('mouseenter.breadcrumb mouseleave.breadcrumb');
            },

            resizeBreadcrumb: function (breadcrumb) {
                var oldbreadcrumbWidth = breadcrumb.data('width');
                var currentBreadcrumbWidth = breadcrumb.width();
                if (oldbreadcrumbWidth < currentBreadcrumbWidth) {
                    this.checkAndExpandBreadCrumb(breadcrumb);
                } else if (oldbreadcrumbWidth > currentBreadcrumbWidth) {
                    this.checkAndCollapseBreadCrumb(breadcrumb);
                }
                breadcrumb.data('width', currentBreadcrumbWidth);
            }
        };

        return this.each(function () {
            //check breadcrumb sizes
            _private.checkBreadcrumb($(this));
            var that = this;
            $(this).bind('window_resize', function() {
                _private.resizeBreadcrumb($(that));
            });
        });
    };
})(jQuery);