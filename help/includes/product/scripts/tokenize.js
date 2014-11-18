(function($) {
  /*
   * Cross-browser setSelectionRange
   *
   * Usage:
   *   $("input[type=text]").setSelectionRange(0, 10).focus();
   */
  $.fn.setSelectionRange = function(start, end) {
    var input = this.get(0);

    if (input.setSelectionRange) {
      input.setSelectionRange(start, end);
    } else if (input.createTextRange) {
      // IE < 9
      var range = input.createTextRange();
      range.collapse(true);
      range.moveEnd("character", end);
      range.moveStart("character", start);
      range.select();
    }
    return this;
  }

  $.fn.cursorPosition = function() {
    var input      = this.get(0),
        dimensions = { start: 0, end: 0 };

    if (input.selectionStart) {
      dimensions.start = input.selectionStart;
      dimensions.end   = input.selectionEnd;
    } else if (document.selection) {
      // IE
      var selection = document.selection.createRange();

      if (selection && selection.parentElement() == input) {
        var range = input.createTextRange();
        range.moveToBookmark(selection.getBookmark());

        var initial_length = range.text.length;

        range.moveStart("character", -input.value.length);

        dimensions.start = range.text.length - initial_length;
        dimensions.end   = range.text.length;
      }
    }
    return dimensions;
  }
})(jQuery);

/*
 * Select elements between two points on an x-axis
 */
(function($) {
  $.fn.between = function(start, end) {
    var matches = $.map(this, function(element) {
      var left  = $(element).offset().left,
          right = left + $(element).width();

      if (right > start && left < end) {
        return element;
      }
    });
    return $(matches);
  }
})(jQuery);

/*
 * Select elements between two points on an x- and y-axis
 */
(function($) {
  $.fn.within = function(start, end) {
    var matches = $.map(this, function(element) {
      var element = $(element),
          left    = element.offset().left,
          right   = left + element.width(),
          top     = element.offset().top,
          bottom  = top + element.outerHeight();

      if (right > start.x && left < end.x &&
          (top <= start.y && bottom >= start.y) &&
          (top <= end.y && bottom >= end.y)) {
        return element.get(0);
      }
    });
    return $(matches);
  }
})(jQuery);

/*
 * tokenize -- turn an text input field into a tokenized widget
 */
(function($) {
  var query_range_helper = function(input) {
    if (input.data("query_range_helper")) {
      return input.data("query_range_helper");
    }

    var helper = $("<span/>").appendTo($("body"));

    helper.css({
      "visibility":     "hidden",
      "position":       "absolute",
      "font-family":    input.css("font-family"),
      "font-size":      input.css("font-size"),
      "font-weight":    input.css("font-weight"),
      "font-style":     input.css("font-style"),
      "letter-spacing": input.css("letter-spacing"),
      "word-spacing":   input.css("word-spacing"),
      "text-indent":    input.css("text-indent"),
      "direction":      input.css("direction"),
      "margin-top":     input.css("margin-top"),
      "margin-right":   input.css("margin-right"),
      "margin-bottom":  input.css("margin-bottom"),
      "margin-left":    input.css("margin-left"),
      "border-top":     input.css("border-top"),
      "border-right":   input.css("border-right"),
      "border-bottom":  input.css("border-bottom"),
      "border-left":    input.css("border-left"),
      "padding-top":    input.css("padding-top"),
      "padding-right":  input.css("padding-right"),
      "padding-bottom": input.css("padding-bottom"),
      "padding-left":   input.css("padding-left")
    });

    input.data("query_range_helper", helper);

    return helper;
  }

  var parse_query = function(query, fields) {
    var result = { query: "", tokens: [] },
        regexp, matchdata;

    $.each(fields, function(i, field) {
      regexp = new RegExp(field + ":(?:\"([^\"]+)\"|'([^']+)'|([^ ]+))", "i");

      while (matchdata = query.match(regexp)) {
        result.tokens.push({
          field: field,
          value: (matchdata[1] || matchdata[2] || matchdata[3])
        });
        query = query.replace(matchdata[0], "");
      }
    });
    result.query = $.trim(query);

    return result;
  }

  var create_token = function(parent, options) {
    var token = $('<span class="token ' + options.field + '"/>');

    token.append($('<span class="label">' + options.label + '</span>'))
         .append($('<span class="remove"><a href="#remove-token">remove</a></span>'))

    token.appendTo(parent).data({
      field: options.field,
      value: options.value
    });

    token.find("a").bind("click", function(event) {
      var form = token.closest("form");
      token.trigger("token.remove");
      form.submit();
      event.preventDefault();
    });
  }

  var tokens_to_query = function(tokens) {
    return $.map(tokens, function(token) {
      return $(token).data("field") + ':' + $(token).data("value");
    }).join(" ");
  }

  $.fn.tokenize = function(options) {
    var input  = this.find("input[type=text]"),
        hidden = $('<input type="hidden" name="' + input.attr("name") + '">'),
        tokens = $('<span class="tokens"/>').prependTo(this),
        parsed_query;

    this.addClass("tokenized");

    // use a hidden input to submit the actual field value
    this.closest("form").append(hidden).bind("submit", function(event) {
      var term = input.val() + " " + tokens_to_query($(event.target).find(".token"));
      term = $.trim(term);

      if (term == "" && options.remove_if_empty) {
        hidden.remove();
      } else {
        hidden.val(term);
      }
    });
    input.removeAttr("name");

    // extract tokens from the query
    parsed_query = parse_query(input.val(), options.fields || []);

    // update the text input
    input.val(parsed_query.query);


    // turn the parsed tokens into elements
    $.each(parsed_query.tokens, function(i, token) {
      create_token(tokens, {
        field: token.field,
        value: token.value,
        label: options.label_function ? options.label_function(token) : token.value 
      });
    });

    // If the user clicks the back button we need to reset the text field
    var reinitialize = function() {
      parsed_query = parse_query(input.val(), options.fields || []);
      input.val(parsed_query.query);
      clearInterval(reinitialize_watcher);
    }

    $(window).unload(function(event) {
      var reinitialize_watcher = setInterval(reinitialize, 1000);
    });

    var self = this;

    this.bind("mousedown", function(event) {
      self.find(".token").removeClass("selected highlighted")
      self.data("drag.start", {x: event.pageX, y: event.pageY});
    });

    // var outline = $("<div/>")
    //       .appendTo($("body"))
    //       .css({
    //         "position": "absolute",
    //         "border": "1px dotted blue",
    //         "height": "100px",
    //         "width": "100px",
    //         "top": "-5000px",
    //         "left": "-5000px"
    //       });

    this.bind("mousemove", function(event) {
      var pos = self.data("drag.start"),
          start, end;

      if (pos) {
        start = {
          x: Math.min(pos.x, event.pageX),
          y: Math.min(pos.y, event.pageY)
        }
        end = {
          x: Math.max(pos.x, event.pageX),
          y: Math.max(pos.y, event.pageY)
        }
        self.find(".token:not(.highlighted)")
            .within(start, end)
            .addClass("highlighted");
      }
    });

    this.bind("mouseup", function(event) {
      var pos = self.data("drag.start"),
          start, end, targets;

      self.data("drag.start", false);

      if (pos) {
        start = {
          x: Math.min(pos.x, event.pageX),
          y: Math.min(pos.y, event.pageY)
        }
        end = {
          x: Math.max(pos.x, event.pageX),
          y: Math.max(pos.y, event.pageY)
        }
        targets = self.find(".token, input")
                      .removeClass("selected highlighted")
                      .within(start, end);

        if (!targets.is("input")) {
          self.find("input").blur();
        }
        if (targets.length > 1 || !targets.is("input")) {
          targets.trigger("token.select", end.x);
        }
      }
    });

    this.bind("token.resize", function(event) {
      self.find("input").width(
        self.width() - self.find(".tokens").width() - 5
      );
    });

    // don't select token text when dragging over it
    this.delegate(".tokens, .token", "selectstart mousedown", function(event) {
      event.preventDefault();
    });

    this.find(".token").bind("token.select", function(event) {
      $(event.target).addClass("selected");
    });

    this.find(".token").bind("token.remove", function(event) {
      $(event.target).remove();
      self.trigger("token.resize");
    });

    this.find("input").bind("token.select", function(event, pageX) {
      var input  = $(event.target),
          helper = query_range_helper(input),
          width  = pageX - input.offset().left,
          end_at = input.val().length,
          index  = 0;

      helper.text(input.val());

      if (width < helper.width()) {
        helper.text("");

        while (index < input.val().length && width > helper.width()) {
          character = input.val().slice(index, index + 1);
          helper.text(helper.text() + character);
          index = index + 1;
        }
        end_at = index;
      }
      input.setSelectionRange(0, end_at).focus();
    });

    this.find("input").bind("keydown", function(event) {
      var input = $(this);

      switch (event.which) {
        // Control-A/Command-A
        case 65:
          if (event.ctrlKey || event.metaKey) {
            self.find(".token").trigger("token.select");
          }
          break;
        // Escape
        case 27:
          self.find(".token").removeClass("selected highlighted");
          input.blur();
          break;
        // Backspace
        case 8:
          if (event.target == input.get(0)) {
            event.stopPropagation();

            input.bind("keyup.token.delete", { "before": input.val() }, function(event) {
              input.unbind("keyup.token.delete");

              if (event.data.before == input.val()) {
                var token = self.find(".token:last");
                if (token) {
                  if (token.hasClass("selected")) {
                    token.trigger("token.remove");
                  } else {
                    token.addClass("selected");
                  }
                }
              } else {
                self.find(".token.selected").trigger("token.remove");
              }
            });
          } 
          break;
      }
    });

    // catch Delete and Backspace when the text field is not focused
    $(document).bind("keydown", function(event) {
      if (event.which == 8 || event.which == 46) {
        self.find(".token.selected").trigger("token.remove");
        if (event.which == 8) {
          event.preventDefault();
        }
      }
    });

    this.trigger("token.resize");

    return this;
  }
})(jQuery);
