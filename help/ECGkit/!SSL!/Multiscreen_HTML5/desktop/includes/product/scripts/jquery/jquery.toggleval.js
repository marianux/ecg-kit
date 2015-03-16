/* -------------------------------------------------- *
 * ToggleVal 3.0
 * Updated: 01/15/2010
 * -------------------------------------------------- *
 * Author: Aaron Kuzemchak
 * URL: http://aaronkuzemchak.com/
 * Copyright: 2008-2010 Aaron Kuzemchak
 * License: MIT License
** -------------------------------------------------- */

(function($) {
	// main plugin function
	$.fn.toggleVal = function(theOptions) {
		// check whether we want real options, or to destroy functionality
		if(!theOptions || typeof theOptions == 'object') {
			theOptions = $.extend({}, $.fn.toggleVal.defaults, theOptions);
		}
		else if(typeof theOptions == 'string' && theOptions.toLowerCase() == 'destroy') {
			var destroy = true;
		}
		
		return this.each(function() {
			// unbind everything if we're destroying, and stop executing the script
			if(destroy) {
				$(this).unbind('focus.toggleval').unbind('blur.toggleval').removeData('defText');
				return false;
			}
			
			// define our variables
			var defText = '';
			
			// let's populate the field, if not default
			switch(theOptions.populateFrom) {
				case 'title':
					if($(this).attr('title')) {
						defText = $(this).attr('title');
						$(this).val(defText);
					}
					break;
				case 'label':
					if($(this).attr('id')) {
						defText = $('label[for="' + $(this).attr('id') + '"]').text();
						$(this).val(defText);
					}
					break;
				case 'custom':
					defText = theOptions.text;
					$(this).val(defText);
					break;
				default:
					defText = $(this).val();
			}
			
			// let's give this field a special class, so we can identify it later
			// also, we'll give it a data attribute, which will help jQuery remember what the default value is
			$(this).addClass('toggleval').data('defText', defText);
			
			// now that fields are populated, let's remove the labels if applicable
			if(theOptions.removeLabels == true && $(this).attr('id')) {
				$('label[for="' + $(this).attr('id') + '"]').remove();
			}
			
			// on to the good stuff... the focus and blur actions
			$(this).bind('focus.toggleval', function() {
				if($(this).val() == $(this).data('defText')) { $(this).val(''); }
				
				// add the focusClass, remove changedClass
				$(this).addClass(theOptions.focusClass);
			}).bind('blur.toggleval', function() {
				if($(this).val() == '' && !theOptions.sticky) { $(this).val($(this).data('defText')); }
				
				// remove focusClass, add changedClass if, well, different
				$(this).removeClass(theOptions.focusClass);
				if($(this).val() != '' && $(this).val() != $(this).data('defText')) { $(this).addClass(theOptions.changedClass); }
					else { $(this).removeClass(theOptions.changedClass); }
			});
		});
	};
	
	// default options
	$.fn.toggleVal.defaults = {
		focusClass: 'tv-focused', // class during focus
		changedClass: 'tv-changed', // class after focus
		populateFrom: 'default', // choose from: default, label, custom, or title
		text: null, // text to use in conjunction with populateFrom: custom
		removeLabels: false, // remove labels associated with the fields
		sticky: false // if true, default text won't reappear
	};
	
	// create custom selectors
	// :toggleval for affected elements
	// :changed for changed elements
	$.extend($.expr[':'], {
		toggleval: function(elem) {
			return $(elem).data('defText') || false;
		},
		changed: function(elem) {
			if($(elem).data('defText') && $(elem).val() != $(elem).data('defText')) {
				return true;
			}
			return false;
		}
	});
})(jQuery);