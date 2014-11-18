//===============================================================================
////==================Search page templates======================================

(function ($) {
    $.setParameter = function(uri, key, value){
        var re = new RegExp("([?&])" + key + "=.*?(&|$)", "i");
        //if the query parameter is already set, simply replace the new value in query parameter.
        if (uri.match(re)) {
            return uri.replace(re, '$1' + key + "=" + value + '$2');
        } else {
            var separator = uri.indexOf('?') !== -1 ? "&" : "?";
            return uri + separator + key + "=" + value;
        }
    };

    $(function() {
        $("#current_year").text(new Date().getFullYear());
    });
})(jQuery);

$.template('searchresults',
        '<div class="searchResult ${$item.getSearchResultClassType}">' +
        '<span class="searchTitle"><a href="${$item.resolveProductHelpDir}/${path}${$item.appendSearchHighlight}">${title}</a></span>' +
        '{{if hassummary == "true"}}<span class="searchSummary"> - ${summary}</span>{{/if}}' +
        '<div class="searchHighlight">{{html $item.getHighlights("...")}}</div>' +
        '<div class="searchBreadcrumb">' +
        '{{each(i, breadcrumb) breadcrumbs}}' +
        '{{each breadcrumb}}' +
            '<a href="${$item.resolveProductHelpDir}/${$value.relativepath}">${$value.label}</a>' +
            '{{if $index < breadcrumb.length-1}} &gt; {{/if}}' +
        '{{/each}}' +
        '{{/each}}</div></div>'
);

$.template('searchsummary', '<div id="resultSummaryTop">${summarytext}</div>');

$.template('footerpages',
    '<div id="search_result_footer">' +
    '{{each(i, summary) searchSummary}}' +
        '<div id="resultSummaryBottom">${summary.summarytext}</div>' +
        '<div id="search_result_pages">' +
        //add the previous label
        '{{if footerdata.selectedpage == 1}}' +
          '<span class="results-page">${$item.getLocalizedString("previous", footerdata.locale)}</span>' +
        '{{else}}' +
          '<a href="${$item.getPageUrl(footerdata.selectedpage-1)}" class="results-page">${$item.getLocalizedString("previous", footerdata.locale)}</a>' +
        '{{/if}}' +
        //if the start range is not 1, add 1 and ...
        '{{if footerdata.startrange != 1}}' +
          '<a href="${$item.getPageUrl(1)}" class="results-page">1</a>' +
          '<span class="results-page">...</span>' +
        '{{/if}}' +
        //add the range
        '{{each  $item.getRangeArray(footerdata.startrange, footerdata.endrange)}}' +
        '{{if $value == footerdata.selectedpage}}' +
           '<span class="results-page-selected">${$value}</span>' +
        '{{else}}' +
           '<a href="${$item.getPageUrl($value)}" class="results-page">${$value}</a>' +
        '{{/if}}' +
        '{{/each}}' +
        //if the end range is less than the total number of pages, add ... and the last page number
        '{{if footerdata.endrange < footerdata.numpages}}' +
         '<span class="results-page">...</span>' +
         '<a href="${$item.getPageUrl(footerdata.numpages)}" class="results-page">${footerdata.numpages}</a>' +
        '{{/if}}'+
        // Add the next label
         '{{if footerdata.selectedpage == footerdata.numpages}}' +
          '<span class="results-page"> ${$item.getLocalizedString("next", footerdata.locale)}</span>' +
        '{{else}}' +
          '<a href="${$item.getPageUrl(footerdata.selectedpage+1)}" class="results-page">${$item.getLocalizedString("next", footerdata.locale)}</a>' +
        '{{/if}}' +
        '</div>' +
    '{{/each}}' +
    '</div>'
);

$.template('facets',
    '{{each(prop, val) $data}}' +
        '<div class="search_refine ${prop}">' +
        '<h3>${val.facettitle}</h3>' +
        '<div class="search_refine_scroll">' +
        '{{if val.hasrefinedfacet=="true"}}' +
        '{{tmpl(val, {facetType: prop}) "refinedfacet"}}' +
        '{{else}}' +
        '{{if val.hasrefinablefacets=="true"}}' +
        '{{tmpl(val, {facetType: prop}) "refinablefacet"}}' +
        '{{/if}}' +
        '{{/if}}' +
        '</div>' +
        '</div>' +
    '{{/each}}'
);

$.template('refinablefacet',
    '<ul>' +
    '{{each(i, result) refinablefacets}}' +
        '<li class="refinable ${result.facetid}_${$item.facetType}">' +
        '<a id="refine-${$item.facetType}-${result.facetid}" href="searchresults.html?${result.faceturlquery}">' +
        '<span class="refine_${$item.facetType}_count">${result.facetcount}</span>' +
        ' ${result.facetname}</a>' +
        '</li>' +
        '{{/each}}' +
    '</ul>'
);

$.template('refinedfacet',
        '<ul>' +
        '<li class="refined ${refinedfacet.facetid}_${$item.facetType}">' +
        '<a href="searchresults.html?${refinedfacet.faceturlquery}">${refinedfacet.facetname}</a>' +
        '{{if refinedfacet.hasrefinedfacet=="true"}}' +
        '{{tmpl(refinedfacet, {facetType: $item.facetType}) "refinedfacet"}}' +
        '{{else}}' +
        '{{if refinedfacet.hasrefinablefacets=="true"}}' +
        '{{tmpl(refinedfacet, {facetType: $item.facetType}) "refinablefacet"}}' +
        '{{/if}}' +
        '{{/if}}' +
        '</li>' +
        '</ul>'
);

$.template('suggestionlist',
    '<div>' +
    '<div><p>${message}</p></div>' +
    '<div><p>${suggestions}</p></div>' +
    '<ul>' +
    '{{each suggestionslist}}' +
        '<li>${$value}</li>' +
    '{{/each}}' +
    '</ul>' +
    '</div>'
);

//===============================================================================

function getSearchResultsHtml(searchResults, searchTerm) {
   return $("<div />").append($.tmpl('searchresults', searchResults, {
       getHighlights: function (separator) {
           return this.data.highlights.join(separator);
       },
       getSearchResultClassType: function() {
           return 'result_type_' + this.data.type;
       },
       resolveProductHelpDir: function() {
           var productHelpDir = this.data.product;
           if (productHelpDir.match(/^file:.*/)) {
               return productHelpDir;
           } else {
               return "../" + productHelpDir;
           }
       },
       appendSearchHighlight: function () {
           if (searchTerm && searchTerm.length > 0) {
               return "?searchHighlight=" + searchTerm;
           } else {
               return "";
           }
       }
    })).html();
}

function getSearchSummaryHtml(searchSummary) {
    return $("<div />").append($.tmpl('searchsummary', searchSummary)).html();
}

function getSearchFooterHtml(searchSummary) {

   return $("<div />").append($.tmpl('footerpages', {

        searchSummary: searchSummary}, {
        getPageUrl: function(pageNumber) {
            return $.setParameter(window.location.href, 'page', pageNumber);
        },

        getRangeArray: function(start, end) {
            var arr = [];
            for (var i = start; i < end; i++) {
                arr.push(i);
            }
            return arr;
        },

        getLocalizedString: function(str, locale) {
            return getLocalizedString(str, locale);
        }
    })).html();
}

function getFacetResultsHtml(facetResults) {
    return $('<div />').append($.tmpl("facets", facetResults)).html();
}

//===============================================================================

function getErrorHtml(error) {
	var html = '';
	html += '<div>';
	html += '<p>';
	html += '<font color="red">';	
	html = html + error;
	html += '</font>';	
	html += '</p>';
	html += '</div>';
	
	return html;
} // end function getErrorHtml

function getSuggestionsListHtml(noResultJson) {
    return $.tmpl("suggestionlist", noResultJson);

}
//===============================================================================
