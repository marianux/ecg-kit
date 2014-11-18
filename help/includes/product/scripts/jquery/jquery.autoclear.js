/**
 * Apply a default value to text fields quickly &amp; easily.
 *
 * The easiest way to use is $('your-selector').autoclear(). All the defaults
 * in the settings object are used. For more advanced cases, and complete
 * reference, @see http://www.mattlunn.me.uk/projects/autoclear
 *
 * @author Matt Lunn
 * @version 6
 * @param  Object / String
 * @return Object jQuery
 * @see http://www.mattlunn.me.uk/projects/autoclear
 * @see README
 */
; (function ($) {
    // Reference to the val function we will be extending
    var _val = jQuery.fn.val;
    var slice = [].slice;

    /**
     * The actual autoclear functionality
     *
     * @author Matt Lunn
     * @param  Object / String
     * @return Object jQuery
     * @see http://www.mattlunn.me.uk/projects/autoclear
     * @see README
     */
    jQuery.fn.autoclear = function (options) {
        var settings = {
            defaultClass: 'default',
            otherClass: 'other',
            defaultValue: '',
            useDefaultOnReset: true,
            clearDefaultOnSubmit: true
        };

        if (arguments.length) {
            switch (typeof options) {
            case "string":
                settings.defaultClass = options;
            break;
            case "object":
                settings = jQuery.extend(settings, options);
            break;
            };
        };

        this.filter('input:text,textarea').each(function () {
            // self caches $(this) for optimization
            var self = $(this);
            // The parent form of the input field.
            var form = self.closest('form');
            // The current value held by the input field
            var currentValue = jQuery.trim(val(self));
            // Holds the value chosen to show as the helper text
            var defaultValue = self.attr('title');

            // Calculate the actual defaultValue. If `title` is empty, we take
            // the current value of the input field unless it's empty, in which
            // case we fall back to the defaultValue in the settings object.
            if (isBlank(defaultValue)) {
                if (currentValue === '') {
                    defaultValue = settings.defaultValue;
                } else {
                    defaultValue = currentValue;
                };
            };

            // Set the default value
            self.data('default.autoclear', defaultValue);			

            // When the form is reset, if useDefaultOnReset is true, we need to
            // set the applied styles of the input field to be default.
            // Otherwise, we need to make sure that the appropiate style is
            // applied to whatever the default value was. self.val() will apply
            // the styles itself.
            form.bind('reset.autoclear', function () {
                var expectedValue;

                if (settings.useDefaultOnReset) {
                    expectedValue = '';
                } else {
                    expectedValue = self.attr('defaultValue');
                };

                self.val(expectedValue);
            });

            // If we have useDefaultOnReset, or there is no defaultValue,
            // set the defaultValue to the helper text. This is so that if a
            // form reset occurs, the intended value is set.
            if (settings.useDefaultOnReset
                || isBlank(jQuery.trim(self.attr('defaultValue')))) {

                self.attr('defaultValue', defaultValue).val(currentValue);
            };

            // If we dont want to submit default values, trigger the focussing
            // of the input on submit; this will clear the default text if it is
            // set.
            if (settings.clearDefaultOnSubmit) {
                form.bind('submit.autoclear', function () {
                    self.trigger('focus.autoclear');
                });
            };
        }).bind({
            'focus.autoclear': function () {
                var self = $(this);

                if (self.val() === '') {
                    val(self, '').trigger('other.autoclear');
                };
            },
            'blur.autoclear': function () {
                var self = $(this);
                var value = jQuery.trim(self.val());

                if (value === '') {
                    self.trigger('default.autoclear');
                } else {
                    self.trigger('other.autoclear');
                };
            },
            'default.autoclear': function () {
                var self = $(this).removeClass(settings.otherClass).
                                   addClass(settings.defaultClass);

                val(self, self.data('default.autoclear'));
            },
            'other.autoclear': function () {
                var self = $(this);

                self.removeClass(settings.defaultClass).
                     addClass(settings.otherClass);
            }
        }).trigger('blur.autoclear');

        return this;
    };

    /**
     * Custom val function. This will still act exactly the same as the native
     * jQuery val function, however if a set operation is carried out with val,
     * the default style/ text will be applied if the value is set to '', and
     * the default style will be removed if the value is changed from '' to a
     * valid value. If a get operation is carried out, and the input is
     * currently showing the default value, '' will be returned instead of the
     * default helper text.
     *
     * @author Matt
     * @see link		http://api.jquery.com/val
     * @param Mixed		Any mixture of valid arguments for native val
     * @return Mixed	Result of native val
     */
    jQuery.fn.val = function () {
        var result = _val.apply(this, arguments);
        var defaultValue;

        // Getter
        if (typeof result === "string") {
            defaultValue = this.first().data    ('default.autoclear');

            if (defaultValue !== undefined && result === defaultValue) {
                result = '';
            };
        } else {
        // Setter

            this.each(function () {
                var self = $(this);
                var defaultValue = self.data('default.autoclear');

                if (defaultValue !== undefined && self.val() === '') {
                    self.trigger('default.autoclear');
                } else {
                    self.trigger('other.autoclear');
                };
            });
        };

        return result;
    };

    /**
     * Helper function which detects whether the provided value is either
     * undefined or an empty string.
     *
     * @author Matt
     * @param 	prop	String	Property to check
     * @return	Boolean	Is it blank?
     */
    function isBlank (prop) {
        return prop === undefined || prop === '';
    };

    /**
     * Helper value which applies the native jQuery val method to the provided
     * jQuery object, rather than applying our custom val method that we replace
     * the native val with. E.g if the value is the helper text, the helper text
     * will be returned, rather than an empty string.
     *
     * @author Matt
     * @param el	Object	jQuery object to apply val to
     * @param arg	Mixed	Optional: Parameter that will be passed to val
     * @return Mixed		Result of applying val to the jQuery Object
     */
    function val (el) {
        return _val.apply(el, slice.call(arguments, 1));
    };

}(jQuery));