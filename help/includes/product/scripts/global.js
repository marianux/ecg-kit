function openWindow ( url, width, height, options, name ) {
  if ( ! width ) width = 640;
  if ( ! height ) height = 420;
  if ( ! options ) options = "scrollbars=yes,menubar=yes,toolbar=yes,location=yes,status=yes,resizable=yes";
  if ( ! name ) name = "outsideSiteWindow";
  var newWin = window.open( url, name, "width=" + width + ",height=" + height + "," + options );
} // end function openWindow