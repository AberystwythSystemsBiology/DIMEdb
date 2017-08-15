$(document).ready(function () {
    $("#autofill_pubmed").click(function () {
        var id = $("input[name=pubmed_id]").val();
        $.getJSON('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id='+id+'&retmode=json', function(summary){
            var author_list = "";
            for(author in summary.result[id].authors){
                author_list += summary.result[id].authors[author].name+', ';
            }


            for (index in summary.result[id].articleids) {
                if (summary.result[id].articleids[index].idtype == "doi") {
                    $("input[name=doi]").val(summary.result[id].articleids[index].value);

                }
            }

            $("input[name=publication_title]").val(summary.result[id].title);
            $("input[name=author_list]").val(author_list);
        });
    });
});
