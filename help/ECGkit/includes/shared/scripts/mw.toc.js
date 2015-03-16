(function($) {
    $.fn.setupToc = function() {
        var _private = {
            fetchTocStrategy: null,
            init: function(toc_container_wrapper) {
                if (this.getTocCookie() === "false") {
                    $('.site_container').removeClass('site_toc_opened').addClass('site_toc_closed');
                }
                var tocIcon = toc_container_wrapper.find('.toc_icon');
                var tocUrl = tocIcon.find('a:first').attr('href');
                tocIcon.data({"tocUrl" : tocUrl});
                tocIcon.data({"globalTocUrl": "../" + tocUrl});
                tocIcon.find('a:first').remove();

                this.addTocHandlers(toc_container_wrapper);
                this.expandToc(toc_container_wrapper);


                tocIcon.click(function(evt) {
                    evt.stopPropagation();
                    _private.toggleToc(toc_container_wrapper);
                });

                this.setTocSize(toc_container_wrapper);


                if ($('#nextgen_b-dc').length > 0) {
                    $(window).scroll(function () {
                        _private.setTocTop(toc_container_wrapper);
                    });
                    _private.setTocTop(toc_container_wrapper);
                }

                var browserInfo = $.browserInfo();
                if (browserInfo.name === 'mozilla' && parseFloat(browserInfo.version.slice(0, 3)) < 3.6) {
                    this.fetchTocStrategy = this.fetchTocForLinux;
                } else {
                    this.fetchTocStrategy = this.fetchTocForAll;
                }

                if ($.fn.nanoScroller) {
                    $('.toc_scroll_pane', toc_container_wrapper).nanoScroller({contentClass: 'toc_container', preventPageScrolling: true});
                    if ($('.current_page').length > 0) {
                        $('.toc_scroll_pane', toc_container_wrapper).nanoScroller({scrollTo:$('.current_page').first()});
                    }


                }

                //Expand the current selected page if possible.
                var currentPageEl = toc_container_wrapper.find('.current_page:first');
                if (currentPageEl && currentPageEl.hasClass('toc_collapsed')) {
                    this.toggleTocListEl(currentPageEl);
                }

            },

            setTocTop: function(toc) {
                var current_top = $(window).scrollTop();
                if (current_top > 0) {
                    if (current_top > 156) {
                        toc.css('top', '10px').css('position', 'fixed');
                    } else {
                        toc.css('top', 'auto').css('position', 'relative');
                    }
                } else {
                    toc.css('top', 'auto').css('position', 'relative')
                }
                toc.data('top', toc.css('top'));
            },

            expandToc: function(tocEl) {
                tocEl.find('li:first').addClass('toc_expanded').removeClass('toc_collapsed');
                $(".current_page:first").parents('li').removeClass('toc_collapsed')
                    .addClass('toc_expanded');
            },

            setTocSize: function(toc_container_wrapper) {
                if (toc_container_wrapper) {
                    var tocHeight;
                    var windowHeight = $(window).height();
                    if ($('#nextgen_b-dc').length > 0) {
                        var webWindowHeight = windowHeight - 196;
                        var contentContainerHeight = $('#content_container').height();
                        tocHeight = Math.min(webWindowHeight, contentContainerHeight);
                    } else {
                        tocHeight = windowHeight - 27;
                    }
                    toc_container_wrapper.css('height', (tocHeight) + 'px').css('max-height',
                        (tocHeight) + 'px');
                    toc_container_wrapper.find('.toc_scroll_pane').css('height',
                        ( tocHeight - 32) + 'px').css('max-height', (tocHeight - 32) + 'px');
                    toc_container_wrapper.find('.toc_container').css('height',
                        ( tocHeight - 37) + 'px').css('max-height', (tocHeight - 37) + 'px');
                }
            },

            fetchToc: function(liEl) {
                if (liEl.children('ul').length > 0) {
                    this.toggleTocListEl(liEl);
                    return;
                }
                this.fetchTocStrategy(liEl);
            },

            shouldCache: function (liEl) {
                if (liEl.data('tocurl') != undefined) {
                    return false;
                }
                return true;
            },

            fetchTocForAll: function (liEl) {
                var filteredId = this.getFilteredId(liEl);

                var tocIcon = $('.toc_icon:first');
                if (tocIcon.data('tocObj') && tocIcon.data('tocObj').length > 0) {
                    var tocObj = tocIcon.data('tocObj');
                    var myContents = $(tocObj.find(filteredId).html()).filter('ul');
                    liEl.append(myContents);
                    this.toggleTocListEl(liEl);
                    return;
                }

                var shouldCache = this.shouldCache(liEl);
                $.get(this.getTocUrl(liEl), {cache: true}, function (data) {
                    var tocObj = $(data);
                    var myContents = $(tocObj.find(filteredId).html()).filter('ul');
                    liEl.append(myContents);
                    if (shouldCache) {
                        tocIcon.data('tocObj', tocObj);
                    }
                    _private.toggleTocListEl(liEl);
                }, 'html');
            },

            fetchTocForLinux: function(liEl) {
                var filteredId = this.getFilteredId(liEl);

                var tocDiv = $('#toc_obj');
                if (tocDiv && tocDiv.length > 0) {
                    var myContents = tocDiv.find(filteredId).children('ul');
                    liEl.append(myContents);
                    this.toggleTocListEl(liEl);
                    return;
                }

                $.get(this.getTocUrl(), {cache: true}, function (data) {
                    var tocDiv = $('<div id="toc_obj" style="display: none;"></div>');
                    tocDiv.get(0).innerHTML = data;
                    $('.toc_pane').append(tocDiv);
                    var myContents = tocDiv.find(filteredId).children('ul');
                    liEl.append(myContents);
                    _private.toggleTocListEl(liEl);
                }, 'html');
            },

            getFilteredId: function (liEl) {
                var id = $(liEl).attr('id');
                return "#" + id.replace(/\./g, '\\.').replace(/\#/g, '\\#');
            },

            addTocHandlers: function(tocEl) {
                tocEl.on('click', 'div.img_wrapper', function (e) {
                    e.stopPropagation();
                    _private.fetchToc($(this).closest('li'));
                });

//                $('li#doc\\.other\\.products').one('click', function (e) {
//                    e.stopImmediatePropagation();
//                    _private.fetchGlobalToc($(this).closest('li'));
//                });

                tocEl.on('click', 'div.text_wrapper', function (e) {
                    e.stopPropagation();
                    if ($(this).find('a').length > 0) {
                        _private.followTocUrl($(this));
                        return false;
                    }
                    _private.fetchToc($(this).closest('li'));
                });

                tocEl.on('mouseenter', "div.text_wrapper",
                    function(e) {
                        if (!tocEl.data('keyboard.navigation')) {
                            var selectedTocListEl = _private.getSelectedTocListEl(tocEl);
                            selectedTocListEl.children('div').children().removeClass('hovering');
                            $(this).addClass('hovering');
                            if ($(this).children('a').length > 0) {
                                _private.correctTocUrl($(this).children('a'));
                                return;
                            }
                            var liParent = $(this).closest('li');
                            if (liParent.hasClass('toc_expanded') ||
                                liParent.hasClass('toc_collapsed')) {
                                $(this).prev('div.img_wrapper').addClass('hovering');
                            }
                        }
                    }).on('mouseleave', 'div.text_wrapper', function(e) {
                        if (!tocEl.data('keyboard.navigation')) {
                            $(this).parent().children().removeClass('hovering');
                        }
                    });


                tocEl.on('mouseenter', 'div.img_wrapper',
                    function(e) {
                        if (!tocEl.data('keyboard.navigation')) {
                            var selectedTocListEl = _private.getSelectedTocListEl(tocEl);
                            selectedTocListEl.children('div').children().removeClass('hovering');

                            var liParent = $(this).closest('li');

                            if (liParent.hasClass('toc_expanded') ||
                                liParent.hasClass('toc_collapsed')) {
                                $(this).addClass('hovering');
                            }
                            if ($(this).next('div.text_wrapper').find('a').length === 0) {
                                $(this).next('div.text_wrapper').addClass('hovering');
                            }
                        }
                    }).on('mouseleave', 'div.img_wrapper', function(e) {
                        if (!tocEl.data('keyboard.navigation')) {
                            $(this).parent().children().removeClass('hovering');
                        }
                    });

                tocEl.on('click', 'div.toc_container', function (e) {
                    e.stopPropagation();
                });

                tocEl.find('div.toc_container_header').click(function (e) {
                    e.stopPropagation();
                });
                tocEl.bind('keydown.toc', function(evt) {
                    var key = evt.keyCode ? evt.keyCode : evt.which;
                    if ((key >= 37 && key <= 40) || key === 13) {
                        _private.handleTocKeyDown(key, tocEl);
                        return false;
                    }
                });
                tocEl.bind('mousemove.toc', function() {
                    tocEl.data('keyboard.navigation', false);
                });
            },

            fetchGlobalToc: function (liEl) {
                $.get(this.getGlobalTocUrl(), {cache: true}, function (data) {
                    var tocObj = $(data);
                    liEl.append(tocObj);
                    _private.toggleTocListEl(liEl);
                }, 'html');
            },

            filterPath: function(string) {
                return string.replace(/(index.html)$/, '')
            },

            followTocUrl: function(textDivElt) {
                var locationPath = this.filterPath(location.pathname);
                if (locationPath.charAt(0) === "/") {
                    locationPath = locationPath.substr(1);
                }
                var anchor = textDivElt.children('a').get(0);
                var anchorPath = anchor.pathname.substring(anchor.pathname.indexOf("#"));
                anchorPath = this.filterPath(anchorPath);
                if (anchorPath.charAt(0) === "/") {
                    anchorPath = anchorPath.substr(1);
                }
                var intrnllnkClicked = false;
                if (locationPath == anchorPath && anchor.hash !== '') {
                    var liEl = textDivElt.closest('li');
                    intrnllnkClicked = true;
                    if (liEl.hasClass('toc_collapsed') || liEl.hasClass('toc_expanded')) {
                        this.fetchToc(liEl);
                        return;
                    }
                }

                //Trigger the toc clicked event.
                $.when(textDivElt.trigger('toc_clicked')).done(function () {
                    if (!textDivElt.children('a').hasClass('corrected_url')) {
                        document.location =
                            this.getTocBaseUrl() + textDivElt.children('a').attr('href');
                    } else {
                        document.location = textDivElt.children('a').attr('href');
                    }
                    if (intrnllnkClicked) {
                        textDivElt.trigger('intrnllnk_clicked');
                    }
                });
            },

            correctTocUrl: function(tocLink) {
                if (!tocLink.hasClass("corrected_url")) {
                    var correctedUrl = this.getTocBaseUrl(tocLink.closest('li')) + tocLink.attr('href');
                    tocLink.attr('href', correctedUrl);
                    tocLink.addClass('corrected_url');
                }
            },

            handleTocKeyDown: function(key, tocEl) {
                tocEl.data('keyboard.navigation', true);
                if (key === 40) {
                    this.handleTocDownKey(tocEl);
                } else if (key === 38) {
                    this.handleTocUpKey(tocEl);
                } else if (key === 37) {
                    this.handleTocLeftKey(tocEl);
                } else if (key === 39) {
                    this.handleTocRightKey(tocEl);
                } else if (key === 13) {
                    this.handleTocEnterKey(tocEl);
                }
            },

            handleTocDownKey: function(tocEl) {
                var selectedTocListEl = this.getSelectedTocListEl(tocEl);
                var nextTocListEl;
                if (selectedTocListEl.length === 0) {
                    nextTocListEl = tocEl.find("li:first")
                } else {
                    if (selectedTocListEl.hasClass('toc_expanded')) {
                        nextTocListEl = selectedTocListEl.children('ul:first').find('li:first');
                    } else {
                        nextTocListEl = this.getNextTocEl(selectedTocListEl);
                    }
                }
                if (nextTocListEl.length > 0) {
                    this.changeTocSelection(nextTocListEl, selectedTocListEl, tocEl);
                }
            },

            handleTocUpKey: function(tocEl) {
                var selectedTocListEl = this.getSelectedTocListEl(tocEl);
                if (selectedTocListEl.length > 0) {
                    var prevTocListEl = this.getPreviousTocEl(selectedTocListEl);
                    if (prevTocListEl.length > 0) {
                        this.changeTocSelection(prevTocListEl, selectedTocListEl, tocEl);
                    }
                }
            },

            handleTocLeftKey: function(tocEl) {
                var selectedTocListEl = this.getSelectedTocListEl(tocEl);
                if (selectedTocListEl.length > 0) {
                    if (selectedTocListEl.hasClass('toc_expanded')) {
                        this.toggleTocListEl(selectedTocListEl);
                    }
                }
            },

            handleTocRightKey: function(tocEl) {
                var selectedTocListEl = this.getSelectedTocListEl(tocEl);
                if (selectedTocListEl.length > 0) {
                    this.fetchToc(selectedTocListEl);
                }
            },

            handleTocEnterKey: function(tocEl) {
                var selectedTocListEl = this.getSelectedTocListEl(tocEl);
                if (selectedTocListEl.length > 0) {
                    var textDiv = selectedTocListEl.children('div').find('div.text_wrapper:first');
                    if (textDiv.find('a').length > 0) {
                        this.followTocUrl(textDiv);
                        return false;
                    }
                    this.toggleTocListEl(selectedTocListEl);
                }
            },

            changeTocSelection: function(newListEl, oldListEl, tocContainer) {
                var textDiv = newListEl.children('div').find('div.text_wrapper');
                var toggleDiv = newListEl.children('div').find('div.img_wrapper');
                textDiv.addClass('hovering');
                if (newListEl.hasClass('toc_expanded') || newListEl.hasClass('toc_collapsed')) {
                    toggleDiv.addClass('hovering');
                }

                newListEl.attr('tabindex', -1);
                newListEl.get(0).focus();
                if (textDiv.children('a').length > 0) {
                    this.correctTocUrl(textDiv.children('a'));
                }
                oldListEl.children('div').find('div').removeClass('hovering');
            },

            getNextTocEl: function(selectedTocEl) {
                var nextTocListEl = [];

                nextTocListEl = selectedTocEl.next("li");
                if (nextTocListEl.length > 0) {
                    return nextTocListEl;
                }

                //If a sibling li element is not found, we need to traverse the tree up.
                //Before traversing the tree up, make sure we do not have any sibling ul elements.
                var siblingUlEl = selectedTocEl.next('ul');
                if (siblingUlEl.length > 0) {
                    nextTocListEl = siblingUlEl.find('li:first');
                }
                if (nextTocListEl.length > 0) {
                    return nextTocListEl;
                }

                // if a sibling li, ul element is not found, traverse the tree up.
                // First check if there are any next sibling 'ul' elements to the parent.
                var parentSiblingUlEl = selectedTocEl.parent('ul').next('ul');
                if (parentSiblingUlEl.length > 0) {
                    nextTocListEl = parentSiblingUlEl.find('li:first');
                }
                if (nextTocListEl.length > 0) {
                    return nextTocListEl;
                }

                //if none of the conditions are met, find the next toc list element recursively using the parent 'li' element as the
                // starting point.
                var parentListEl = selectedTocEl.parents('li:first');
                if (parentListEl.length > 0) {
                    return this.getNextTocEl(parentListEl);
                }
                return nextTocListEl;
            },

            getPreviousTocEl: function(selectedTocEl) {
                var prevTocListEl = [];

                // It is important that these checks are performed in this order.
                //First, check if there exists a sibling 'li' element.
                prevTocListEl = selectedTocEl.prev('li');
                if (prevTocListEl.length > 0) {
                    //If a previous sibling does exist, make sure we check if the li element is expanded.
                    //If the previous sibling is expanded, get the last open/visible toc list element in that subtree.
                    prevTocListEl = this.getLastOpenTocListEl(prevTocListEl);
                    return prevTocListEl;
                }

                //if a sibling li element is not found, we need to traverse the tree upwards/backwards
                // before traversing the tree up, make sure we do not have any sibling ul elments;
                var siblingUlEl = selectedTocEl.prev('ul');
                if (siblingUlEl.length > 0) {
                    prevTocListEl = siblingUlEl.children('li:last');
                    prevTocListEl = this.getLastOpenTocListEl(prevTocListEl);
                }
                if (prevTocListEl.length > 0) {
                    return prevTocListEl;
                }

                // if a sibling li, ul element is not found, traverse the tree up.
                // First check if there are any previous sibling 'ul' elements to the parent.
                var parentSiblingUlEl = selectedTocEl.parent('ul').prev('ul');
                if (parentSiblingUlEl.length > 0) {
                    prevTocListEl = parentSiblingUlEl.children('li:last');
                    prevTocListEl = this.getLastOpenTocListEl(prevTocListEl);
                }
                if (prevTocListEl.length > 0) {
                    return prevTocListEl;
                }

                // if the selected toc list element is the first in the list, find the parent toc list element.
                // Do not try to check if the toc list element is expanded, as that would lead to an infinite loop.
                var selectedTocElIndex = selectedTocEl.index();
                if (selectedTocElIndex === 0) {
                    prevTocListEl = selectedTocEl.parents("li:first");
                }

                if (prevTocListEl.length > 0) {
                    return prevTocListEl;
                }

                //if none of the conditions are met, find the previous toc list element recursively using the parent 'li' element as the
                // starting point.
                var parentListEl = selectedTocEl.parents('li:first').prev("li");
                if (parentListEl.length > 0) {
                    return this.getPreviousTocEl(parentListEl);
                }
                return prevTocListEl;
            },

            getLastOpenTocListEl: function(tocListEl) {
                var isExpanded = true;
                while (isExpanded) {
                    if (tocListEl.hasClass('toc_expanded')) {
                        tocListEl = tocListEl.children('ul:last').children('li:last');
                    } else {
                        isExpanded = false;
                    }
                }
                return tocListEl;
            },

            getSelectedTocListEl: function(tocEl) {
                return tocEl.find('.hovering:first').closest('li');
            },

            toggleToc: function(tocEl) {
                if (tocEl.find('.toc_container').is(":visible")) {
                    $('.site_container').removeClass('site_toc_opened').addClass('site_toc_closed');
                    this.setTocCookie('false');
                    tocEl.find('.toc_scroll_pane').nanoScroller({stop: true});
                } else {
                    $('.site_container').removeClass('site_toc_closed').addClass('site_toc_opened');
                    this.setTocCookie('true');
                    tocEl.find('.toc_scroll_pane').nanoScroller({start: true});
                }
                tocEl.trigger('toc_resize');
            },

            toggleTocListEl: function(liEl) {
                if (liEl.hasClass('toc_expanded')) {
                    liEl.children('ul').slideUp('fast', function () {
                        liEl.removeClass('toc_expanded').addClass('toc_collapsed');
                        $('.toc_scroll_pane').nanoScroller({reset: true});
                    });
                } else if (liEl.hasClass('toc_collapsed')) {
                    liEl.children('ul').slideDown('fast', function () {
                        liEl.removeClass('toc_collapsed').addClass('toc_expanded');
                        $('.toc_scroll_pane').nanoScroller({reset: true});
                    });
                }
            },

            getTocUrl: function (liEl) {
                if (liEl && liEl.data('tocurl') != undefined) {
                    if (liEl.parents("#doc\\.other\\.products").length > 0) {
                        return this.getGlobalTocBaseUrl() + liEl.data('tocurl');
                    } else {
                        return liEl.data('tocurl');
                    }
                }
                return $('.toc_icon:first').data('tocUrl');
            },

            getGlobalTocBaseUrl: function () {
                var globalTocUrl = $('.toc_icon:first').data('globalTocUrl');
                return globalTocUrl.substring(0, globalTocUrl.indexOf("doccentertoc"));
            },

            //Get the correct base URL for the TOC links, based on the current page.
            getTocBaseUrl: function(liEl) {
                var tocUrl = this.getTocUrl(liEl);
                return tocUrl.substring(tocUrl, tocUrl.indexOf("doccentertoc"));
            },

            getTocCookie: function() {
                var tocFlag = "MW_toc_docked";
                if (this.isLocalStorageEnabled()) {
                    return window.localStorage.getItem(tocFlag);
                } else {
                    var tocCookie = tocFlag + "=";
                var cookies = document.cookie.split(';');
                for (var i = 0; i < cookies.length; i++) {
                    var cookie = $.trim(cookies[i]);
                        if (cookie.indexOf(tocCookie) === 0) {
                            return cookie.substring(tocCookie.length, cookie.length);
                        }
                    }
                }
                return null;

            },

            setTocCookie: function(value) {
                var tocFlag = "MW_toc_docked";
                if (this.isLocalStorageEnabled()) {
                    window.localStorage.setItem(tocFlag, value);
                } else {
                var date = new Date();
                date.setTime(date.getTime() + (7 * 24 * 60 * 60 * 1000));
                var expiresDate = date.toGMTString();
                    document.cookie = tocFlag + "=" + value
                    + "; expires=" + expiresDate
                    + "; path=/";
                }
            },

            isLocalStorageEnabled: function () {
                return !!window.localStorage;
            }
        };

        return this.each(function () {
            _private.init($(this));
            var that = this;
            $(this).bind('window_resize', function() {
                _private.setTocSize($(that));
            });
        });
    };
})(jQuery);
