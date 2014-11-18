//===============================================================================
////==================HSP link templates=========================================

$.template('installedlink',
        '<li class="installed">' +
            '<a href="${helplocation}">${displayname}</a>' +
        '</li>'
);

$.template('weblink',
        '<div id="weblink">' +
            '<a href="${helplocation}">${displayname}</a>' +
        '</div>');

//===============================================================================

function getInstalledLinkHtml(installedlinks) {
    var hspListHtml =  $('<div/>').append($.tmpl("installedlink", installedlinks)).html();
    return '<ul>' + hspListHtml +'</ul>';
}

function getWebLinkHtml(weblink) {
    return $.tmpl('weblink', weblink).html();
}

//===============================================================================
