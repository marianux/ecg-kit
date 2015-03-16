// Copyright 2013 MathWorks, Inc.

function StringMap() {
    this.strings = {
        en: {
            click_to_collapse:'Click to Collapse',
            click_to_expand:'Click to Expand',
            back_to_top_of_page:'Back to Top of Page',
            back_to_top_of_section:'Back to Top of Section',
            yes:'Yes',
            no:'No',
            was_this_topic_helpful:'Was this topic helpful?',
            expand_all:'expand all',
            expand_all_in_page:'expand all in page',
            collapse_all:'collapse all',
            collapse_all_in_page:'collapse all in page',
            play:'Play',
            stop:'Stop',
            search_suggestions:'Search Suggestions',
            next: 'Next',
            previous: 'Previous'
        },
        ja_JP: {
            click_to_collapse:'クリックして折りたたむ',
            click_to_expand:'クリックして展開する',
            back_to_top_of_page:'ページのトップへ',
            back_to_top_of_section:'セクションのトップへ',
            yes:'はい',
            no:'いいえ',
            was_this_topic_helpful:'<span style=\"font-size:1.2em\">この情報は役に立ちましたか?</span>',
            expand_all:'すべて展開する',
            expand_all_in_page:'ページ内をすべて展開する',
            collapse_all:'すべて折りたたむ',
            collapse_all_in_page:'ページ内をすべて折りたたむ',
            play:'再生',
            stop:'停止',
            search_suggestions:'検索文字列の候補',
            next: '次へ',
            previous: '前へ'
        },
        ko_KR: {
            click_to_collapse:'축소하려면 클릭하십시오',
            click_to_expand:'확장하려면 클릭하십시오',
            back_to_top_of_page:'페이지 맨 위로 돌아가기',
            back_to_top_of_section:'섹션 맨 위로 돌아가기',
            yes:'예',
            no:'아니요',
            was_this_topic_helpful:'이 항목이 도움이 되셨습니까?',
            expand_all:'모두 확장',
            expand_all_in_page:'페이지 내 모두 확장',
            collapse_all:'모두 축소',
            collapse_all_in_page:'페이지 내 모두 축소',
            play:'재생',
            stop:'중지',
            search_suggestions:'검색 제안',
            next: '다음',
            previous: '이전'
        },
        zh_CN: {
            click_to_collapse:'点击以折叠',
            click_to_expand:'点击以展开',
            back_to_top_of_page:'返回页首',
            back_to_top_of_section:'返回节首',
            yes:'是',
            no:'否',
            was_this_topic_helpful:'<span style=\"font-size:1.2em\">本主题对您是否有帮助？</span>',
            expand_all:'全部展开',
            expand_all_in_page:'全页展开',
            collapse_all:'全部折叠',
            collapse_all_in_page:'全页折叠',
            play:'播放',
            stop:'停止',
            search_suggestions:'﻿搜索建议',
            next: '下一页',
            previous: '上一页'
        }
    };

    this.strings['ja'] = this.strings['ja_JP'];
    this.strings['ko'] = this.strings['ko_KR'];
    this.strings['zh'] = this.strings['zh_CN'];

    this.getLocalizedString = function(lang, str) {
        return this.strings[lang][str];
    };
}

function getPageLanguage() {
    return $('.toc_header').attr('lang');
}

function getLocalizedString(str, locale) {
    var lang = locale ? locale : getPageLanguage() || 'en';
    var sMap = new StringMap();
    return sMap.getLocalizedString(lang, str);
}