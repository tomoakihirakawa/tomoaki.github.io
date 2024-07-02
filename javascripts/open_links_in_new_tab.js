document.addEventListener("DOMContentLoaded", function() {
    // 画像リンクに対してtarget="_blank"を追加
    var imageLinks = document.querySelectorAll('a[href$=".png"], a[href$=".jpg"], a[href$=".jpeg"], a[href$=".gif"]');
    imageLinks.forEach(function(link) {
        link.setAttribute('target', '_blank');
    });

    // 外部リンクに対してtarget="_blank"を追加
    var allLinks = document.querySelectorAll('a');
    allLinks.forEach(function(link) {
        if (link.hostname !== window.location.hostname) {
            link.setAttribute('target', '_blank');
        }
    });
});  
