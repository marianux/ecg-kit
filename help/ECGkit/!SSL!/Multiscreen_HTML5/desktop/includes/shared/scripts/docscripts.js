// Toggle a division as expanded or collapsed.
// Also toggle the arrow icon.
// Refer to the division and image by their IDs.
//
// "Collapsed" material is hidden using the
// display property in CSS.

// Used by adaptproduct function (see below)
// to support adaptive doc in the Windows
// version of the Help Browser.
var adaptiveIds = new Array();
function toggleexpander(blockid, arrowid) {
    arrow = document.getElementById(arrowid);
    block = document.getElementById(blockid);
    if (block.style.display == "none") {
        // Currently collapsed, so expand it.
        block.style.display = "block";
        if (arrow.src.indexOf("/") >= 0) { 
           arrow.src = arrow.src.substring(0,arrow.src.lastIndexOf("/")+1) + "arrow_down.gif";
        } else {
           arrow.src = "arrow_down.gif";
        }
        arrow.title = getLocalizedString('click_to_collapse');
    }
    else {
        // Currently expanded, so collapse it.
        block.style.display = "none";
        if (arrow.src.indexOf("/") >= 0) { 
           arrow.src = arrow.src.substring(0,arrow.src.lastIndexOf("/")+1) + "arrow_right.gif";
        } else {
           arrow.src = "arrow_right.gif";
        }
        arrow.title = getLocalizedString('click_to_expand');
    }
    return false; // Make browser ignore href.
}

// ===================================================
// Create and uniquely name two levels of upward navigation buttons
// for Functions -- By Category pages

var top_button_count = 0;
var current_section_id = 0;

function addTopOfPageButtons() {

    top_button_count = top_button_count + 1;

    var top_of_page_buttons =

        "<a class=\"pagenavimglink\" href=\"#top_of_page\" onMouseOver=\"document.images.uparrow" +
            top_button_count +
            ".src=\'doc_to_top_down.gif\'\;\" onMouseOut=\"document.images.uparrow" +
            top_button_count +
            ".src=\'doc_to_top_up.gif\'\;\"><img style=\"margin-top:0;margin-bottom:0px;padding-top:0;padding-bottom:0\" border=0 src=\"doc_to_top_up.gif\"  alt=\"" +
            getLocalizedString('back_to_top_of_page') +
            "\" title=\"" +
            getLocalizedString('back_to_top_of_page') +
            "\" name=\"uparrow" +
            top_button_count +
            "\">\&nbsp\;</a>";

    document.write(top_of_page_buttons);
}


function updateSectionId(id) {
    current_section_id = id;
}


function addTopOfSectionButtons() {

    top_button_count = top_button_count + 1;

    var top_of_page_buttons =

        "<a class=\"pagenavimglink\" href=" +
            "\"#" + current_section_id + "\"" +
            " onMouseOver=\"document.images.uparrow" +
            top_button_count +
            ".src=\'doc_to_section_down.gif\'\;\" onMouseOut=\"document.images.uparrow" +
            top_button_count +
            ".src=\'doc_to_section_up.gif\'\;\"><img style=\"margin-top:0;margin-bottom:0px;padding-top:0;padding-bottom:0\" border=0 src=\"doc_to_section_up.gif\"  alt=\"" +
            getLocalizedString('back_to_top_of_section') +
            "\" title=\"" +
            getLocalizedString('back_to_top_of_section') +
            "\" name=\"uparrow" +
            top_button_count +
            "\">\&nbsp\;</a>";

    document.write(top_of_page_buttons);
}

// ===================================================
// Create and write to the document stream HTML for 
// the link to the Doc Feedback Survey site.
//
// Doing this through a JavaScript function is necessary
// to work around the an issue with pages that are found
// through the search facility of the help browser--
//
// When found as the result of a search, 
// the document that is displayed in the Help browser
// is actually a temporary document with a trivial URL
// such as "text://5", not an actual page location.
//
// But the Help browser inserts a <BASE> element at the beginning
// of each such temporary page, and the <BASE> element stores the
// actual location. 
//
// So this function tests the URL of the document for the expression "text://"
// and if that expression is found, attempts to use the URL stored in
// the <BASE> element.

function writeDocFeedbackSurveyLink() {
    var queryexpression = document.location.href;

    if (queryexpression.search(/text:\/\//) != -1) {
        var baseelement = document.getElementsByTagName("BASE")[0];
        queryexpression = baseelement.href;
    }
    survey_url_yes = "http://www.customersat3.com/TakeSurvey.asp?si=YU2FDmNEifg%3D&SF=" + queryexpression + "-YES";
    survey_url_no = "http://www.customersat3.com/TakeSurvey.asp?si=YU2FDmNEifg%3D&SF=" + queryexpression + "-NO";

    code = '<div style="padding-right:10px" class="feedbackblock">' + '<strong>' + getLocalizedString('was_this_topic_helpful') + '</strong> <input type="button" value="' + getLocalizedString('yes') + '" onClick="openWindow(\'' + survey_url_yes + '\',850,680, \'scrollbars=yes,resizable=yes\'); return false;"/>' + '&nbsp;&nbsp;' + '<input type="button" value="' + getLocalizedString('no') + '" onClick="openWindow(\'' + survey_url_no + '\',850,680, \'scrollbars=yes,resizable=yes\'); return false;"/>' + '</div>';
    document.write(code);
}


// Utility function replacing openWindow function used by the web-site survey link code.
// In the help browser, the original code would create a blank window before loading the URL into the system browser.
function openWindow(url, width, height, options, name) {
    // ignore the arguments, except url
    document.location = url;
} // end function openWindow


// Utility function for linking to feedback survey, as of R2012b.
function openFeedbackWindow(url) {
    window.open(url,"_blank");
} // end function openFeedbackWindow



// ===================================================
// Workaround for G801125.
// This global object check tests for IE8 or lower.
if (document.all && !document.getElementsByClassName) {
    document.createElement("section");
}



// ===================================================
// Function reference pages

$(window).load(function () {
    // Perform breadcrumb check in window load, since all the images in the breadcrumb
    // need to be loaded for correct width calculations.
    getBreadcrumb().mwBreadcrumb();

    //delay the expanding to enable the page to load completely.
    setTimeout(function() {
        expandCollapsedContent();
    }, 0);
});

$(document).ready(function () {
    // this function is derived from an earlier version of jquery
    // and we are using it here for backwards compatability.
    $.browserInfo = function() {
        var ua = navigator.userAgent.toLowerCase();
        var match = /(chrome)[ \/]([\w.]+)/.exec(ua) ||
            /(webkit)[ \/]([\w.]+)/.exec(ua) ||
            /(opera)(?:.*version|)[ \/]([\w.]+)/.exec(ua) ||
            /(msie) ([\w.]+)/.exec(ua) ||
            ua.indexOf("compatible") < 0 && /(mozilla)(?:.*? rv:([\w.]+)|)/.exec(ua) ||
            [];

        var browserName = match ? match [1] : navigator.appName;
        var browserVersion = match ? match[2] : navigator.appVersion;
        return {
            name: browserName,
            version: browserVersion
        }
    };

    if (getParameterByName("browser") === "F1help") {
        $('.site_container:first').removeClass('site_toc_closed').addClass('navigation_off');
        $("#doc_center_content").on('mouseenter', 'a', function() {
            if ($(this).attr('href').indexOf("matlab:") === 0) {
                return;
            }
            if ($(this).hasClass('corrected_url')) {
                return;
            }
            $(this).addClass('corrected_url');
            $(this).attr('href', function (i, h) {
                if (h === undefined) {
                    return "";
                }
                var srcUrl, hash, hashIndex;
                if (h.indexOf('#') > 0) {
                    hashIndex = h.indexOf('#');
                    hash = h.substring(hashIndex, h.length);
                } else {
                    hash = '';
                    hashIndex = h.length;
                }

                srcUrl = h.substring(0, hashIndex);
                return srcUrl + (srcUrl.indexOf('?') != -1 ? "&browser=F1help" : "?browser=F1help") + hash;
            });
        });
    } else {
        $('div.toc_container_wrapper').setupToc();
    }

    //Perform JS code which has any user visible impact first.
    //Check image sizes. Do not scale animated images or any images with hotspots.


    $('#doc_center_content img:not(".animated-image, [usemap]"), #content_container2 img:not(".animated-image, [usemap]")').scaleImage();
    $('#doc_center_content img.animated-image, #content_container2 img.animated-image').animateImage();

    addSmoothScroll();

    $('#content_container .expandAllLink').click(function(e) {
        e.stopPropagation();
        var link = $(this);
        if (link.data('allexpanded')) {
            doCollapse(link);
            setExpandAllLinkState(link, "collapsed");
        } else {
            doExpand(link);
            setExpandAllLinkState(link, "expanded");
        }
    });

    $("#content_container").delegate(".expand", "click", function(e) {
        e.stopPropagation();
        doToggle($(this));
        return false;
    });

    $('#expandAllPage').click(function () {
        if ($(this).data('allexpanded')) {
            doCollapseAllInPage($('#content_container .expand'));
            $(this).data('allexpanded', false);
            //todo: localize
            $(this).html(getLocalizedString('expand_all_in_page'));
        } else {
            doExpandAllInPage($('#content_container .expand'));
            $(this).data('allexpanded', true);
            //todo: localize
            $(this).html(getLocalizedString('collapse_all_in_page'));
        }
    });
    applySearchHighlight();

    $(window).bind('toc_resize', function() {
        $('#search_crumb_container').width($('#content_container').width());
        getBreadcrumb().trigger('window_resize');
    });

    $(window).bind('resize', function(e) {
        $('#search_crumb_container').width($('#content_container').width());
        getBreadcrumb().trigger('window_resize');
        $('.toc_container_wrapper').trigger('window_resize');
    });

    $(window).bind('content_resize', function (e) {
        $('.toc_container_wrapper').trigger('window_resize');
    });

      $(window).bind('intrnllnk_clicked', function (e) {
          expandCollapsedContent();
      });

    $('#search_crumb_container').width($('#content_container').width());
    getBreadcrumb().trigger('window_resize');
});

function setExpandAllLinkState(link, state) {
    if (state === 'expanded') {
        link.data('allexpanded', true);
        link.html(getLocalizedString('collapse_all'));
    } else if (state === 'collapsed') {
        link.data('allexpanded', false);
        link.html(getLocalizedString('expand_all'));
    }
}

function applySearchHighlight() {
    if ($.fn.highlight) {
    var highlighterCSSClass = ['highlight_01', 'highlight_02' ,'highlight_03', 'highlight_04', 'highlight_05'];
        var searchHighlightTerm = getParameterByName('searchHighlight');
        if (searchHighlightTerm.length > 0) {
            var searchHighlightArray = getParameterByName('searchHighlight').match(/"[^"]+"|\S+/g);
            $.each(searchHighlightArray, function (index, value) {
                var searchTerm = value.replace(/^"|"$/g, '');
                var cssClass = highlighterCSSClass[index % highlighterCSSClass.length];
                $("#doc_center_content").highlight(searchTerm,
                    {className:cssClass, wordsOnly: true});
                var elements = $("#doc_center_content").find("." + cssClass);

                setTimeout(function() {
                    expandHighlightedElements(elements);
                }, 50);
            });

            $(document).keyup(function (e) {
                if (e.which === 27) {
                    var classArray = $.map(highlighterCSSClass, function(value) {
                        return "." + value;
                    });
                    var highlightedEl = $(classArray.join(","));
                    $.each(highlightedEl, function() {
                        $(this).removeClass();
                    });
                }
            });
        }
    }
}

function expandHighlightedElements(elements) {
    $.each(elements, function() {
        if (!$(this).is(":visible")) {
            var collapsedParent = $(this).closest('.collapse').prev();
            if (collapsedParent.length > 0) {
                prepareEltForExpansion(collapsedParent, true);
                if (collapsedParent.hasClass('expand') &&
                    !collapsedParent.hasClass('expanded')) {
                    doToggle(collapsedParent, true);
                }
            }
        }
    });
}

function getParameterByName(name) {
    name = name.replace(/[\[]/, "\\\[").replace(/[\]]/, "\\\]");
    var regexS = "[\\?&]" + name + "=([^&#]*)";
    var regex = new RegExp(regexS, 'g');
    var results;
    var value = "";
    while (true) {
        results = regex.exec(window.location.href);
        if (results == null) {
            break;
        }
        value = decodeURIComponent(results[1].replace(/\+/g, " "));
    }
    return value;
}


//helper method to fetch the breadcrumb.
function getBreadcrumb() {
    var breadcrumb;
    if ($("#breadcrumbs").length != 0) {
        breadcrumb = $("#breadcrumbs");
    } else {
        breadcrumb = $(".breadcrumbs:first");
    }
    return breadcrumb;
}

function expandCollapsedContent() {
    if (location.hash.length > 0) {

        var target = getInternalLinkTarget(location.hash);
        if (target.length > 0) {
            var nextSibling = getNextSiblingForAnchorTarget(target);
            prepareEltForExpansion(nextSibling);

            //scroll to the target first.
            var scrollParameter = getScrollParameter();
            var scrollTop = target.offset().top - getScrollTopAdjustment();
            $(scrollParameter).scrollTop(scrollTop);

            if (nextSibling.hasClass('expand') && !nextSibling.hasClass('expanded')) {
                doToggle(nextSibling);
                nextSibling.addClass('anchor_hinting');
                setTimeout(function () {
                    nextSibling.removeClass('anchor_hinting');
                }, 5000);
            }
        }
    }
}

function prepareEltForExpansion(elt, noAnimation) {
    if (elt.hasClass('expand')) {
        doExpandNestedParent(elt, noAnimation);
    } else if (!elt.is(":visible")) {
        doExpandParent(elt, noAnimation);
    }
}

function getNextSiblingForAnchorTarget(target) {
    var nextSibling;
    if (target.is('a')) {
        nextSibling = target.next();
    } else {
        //This needs to be cleaned up.
        //We need to make sure all anchor targets are anchor tags themselves.
        //
        if (target.children().length > 0) {
            nextSibling = target.children(":first");
        } else {
            nextSibling = target;
        }
    }
    return nextSibling;
}

function addSmoothScroll() {
    $(".intrnllnk").each(function() {
        var hash = this.hash;
        var target = getInternalLinkTarget(hash);
        if (target.length > 0) {
            $(this).click(function (evt) {
                evt.preventDefault();
                var nextSibling = getNextSiblingForAnchorTarget(target);
                prepareEltForExpansion(nextSibling);

                var scrollParameter = getScrollParameter();
                var scrollTop = target.offset().top - getScrollTopAdjustment();

                $(scrollParameter).animate({scrollTop:scrollTop}, 700, function () {
                    nextSibling.addClass('anchor_hinting');
                    setTimeout(function () {
                        nextSibling.removeClass('anchor_hinting');
                    }, 5000);
                    if (nextSibling.hasClass('expand') && !nextSibling.hasClass('expanded')) {
                        doToggle(nextSibling);
                    }
                });
                location.hash = hash;
            })
        }
    });
}

function getInternalLinkTarget(hash) {
    //search for anchor with given hash as "name" atrribute value;
    var target = [];

    //Remove the first '#' character from the name attribute. Escape any special character from the name/id.
    var escapedHash = hash.substring(1).replace(/([;&,.+*~':"!^#$%@\[\]\(\)=>\|])/g, '\\$1');

    target = $('[name=' + escapedHash + ']');
    //If no target is found, try to find an element whose id has value = hash.
    if (target.length === 0) {
        target = $(hash);
    }
    return target;
}

function findExpandableContent(elt) {
    if (!elt.hasClass('expandableContent')) {
        elt = elt.closest('.expandableContent');
    }
    return elt;
}

function triggerContentResize() {
    $(window).trigger('content_resize');
}

function doExpand(elt, noAnimation) {
    var expandable = findExpandableContent(elt);
    var browserInfo = $.browserInfo();
    if ((browserInfo.name === 'msie' && browserInfo.version <= 8) || noAnimation) {
        expandable.find('.collapse').show();
        expandable.find('.expand').addClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
        return;
    }

    expandable.find('.collapse').slideDown(function () {
        expandable.find('.expand').addClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
    });
}

function doCollapse(elt, noAnimation) {
    var expandable = findExpandableContent(elt);
    //Before collapsing, check if the collapse child has any expandableContent children.
    //If it does, those divs need to be collapsed and not the parent.
    var collapsedChild = expandable.children('.collapse');
    if (collapsedChild.children('.expandableContent').length > 0) {
        expandable = expandable.children('.collapse').children('.expandableContent');
    }

    var browserInfo = $.browserInfo();
    if ((browserInfo.name === 'msie' && browserInfo.version <= 8) || noAnimation) {
        expandable.find('.collapse').hide();
        expandable.find('.expand').removeClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
        return;
    }

    expandable.find('.collapse').slideUp(function () {
        expandable.find('.expand').removeClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
    });
}

function doToggle(elt, noAnimation) {
    var expandable = findExpandableContent(elt);
    var browserInfo = $.browserInfo();
    if ((browserInfo.name === 'msie' && browserInfo.version <= 8) || noAnimation) {
        expandable.children('.collapse').toggle();
        expandable.children('.expand').toggleClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
        return;
    }

    expandable.children('.collapse').slideToggle(function () {
        expandable.children('.expand').toggleClass('expanded');
        triggerContentResize();
        checkExpandAllLinkState(elt);
    });
}

function checkExpandAllLinkState(elt) {
    //Check if the expandable elt is nested within another expandable elt.
    var expandableParent = elt.parents('.expandableContent:eq(1)');

    // If element is not nested, or there is not expand all link, return
    if (expandableParent.length === 0 || expandableParent.find('.expandAllLink').length === 0) {
        return;
    }
    var expandAllLink = expandableParent.find('.expandAllLink:first');
    var expandableChildren = expandableParent.find('.expandableContent');


    if (elt.hasClass('expanded')) {
        var allChildrenExpanded = true;
        expandableChildren.each(function() {
            var expand = $(this).find('.expand');
            if (!expand.hasClass("expanded")) {
                allChildrenExpanded = false;
            }
        });
        if (allChildrenExpanded) {
            setExpandAllLinkState(expandAllLink, "expanded");
        }
    } else {
        var allChildrenCollapsed = true;
        expandableChildren.each(function() {
            var expand = $(this).find('.expand');
            if (expand.hasClass("expanded")) {
                allChildrenCollapsed = false;
            }
        });

        if (allChildrenCollapsed) {
            setExpandAllLinkState(expandAllLink, "collapsed");
        }
    }
}

function doExpandNestedParent(elt, noAnimation) {
    var expandable = findExpandableContent(elt);
    var expandableParent = expandable.parent().siblings('.expand');
    if (expandableParent.length > 0) {
        if (!expandableParent.hasClass('expanded')) {
            doToggle(expandableParent, noAnimation);
        }
    }
}

// This method is used to support the old HTML template for the ref pages. When all the ref pages have been
// converted to the new template, consolidate this method with the doExpandNestedParent method above.
function doExpandParent(elt, noAnimation) {
    var expandable = elt.closest('.collapse');
    var expandableParent = expandable.siblings('.expand');
    if (expandableParent.length > 0) {
        if (!expandableParent.hasClass('expanded')) {
            doToggle(expandableParent, noAnimation);
        }
    }
}

function doExpandAllInPage(elt) {
    var expandable = findExpandableContent(elt);
    expandable.find('.collapse').show();
    expandable.find('.expand').addClass('expanded');
    triggerContentResize();
    $.each($('#content_container .expandAllLink'), function(i, link) {
        setExpandAllLinkState($(link), "expanded");
    });
}

function doCollapseAllInPage(elt) {
    var expandable = findExpandableContent(elt);
    expandable.find('.collapse').hide();
    expandable.find('.expand').removeClass('expanded');
    triggerContentResize();
    $.each($('#content_container .expandAllLink'), function(i, link) {
        setExpandAllLinkState($(link), "collapsed");
    });

}

function getScrollParameter() {
    //On IE and FF, the slow scroll parameter is the HTML dom element. On webkit, it is the body.
    var browserInfo = $.browserInfo();
    var scrollParameter = (browserInfo.name === 'msie' || browserInfo.name === 'mozilla') ?
        'html' : 'body';
    return scrollParameter;
}

function getScrollTopAdjustment() {
    var scrollTop = 0;
    if ($('#search_crumb_container').css('position') === 'fixed') {
        scrollTop = 103;
    } else {
        scrollTop = 10;
    }
    return scrollTop;
}

// Copyright 2002-2013 The MathWorks, Inc.
