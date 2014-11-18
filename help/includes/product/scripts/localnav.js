$(document).ready(function(){
  $('li.expandable > a').click(function() {
    $(this).parent().toggleClass("collapsed").toggleClass("expanded"); return false;
  });	
});