function getBaseURL () {
   return location.protocol + "//" + location.hostname +
      (location.port && ":" + location.port) + "/";
}


function clear_card() {
    $("#card_metabolite_name").empty();
    $("#molecular_formula").empty();
    $("#neutral_mass").empty();
    $("#positive_adduct").empty();
    $("#negative_adduct").empty();
    $("#origins").empty();
    $("#biofluids").empty();
    $("#tissues").empty();
}

$(document).on('click', '.panel-heading span.clickable', function(e){
    var $this = $(this);
	if(!$this.hasClass('panel-collapsed')) {
		$this.parents('.panel').find('.panel-body').slideDown();
		$this.addClass('panel-collapsed');
		$this.find('i').removeClass('glyphicon-chevron-up').addClass('glyphicon-chevron-down');
	} else {
		$this.parents('.panel').find('.panel-body').slideUp();
		$this.removeClass('panel-collapsed');
		$this.find('i').removeClass('glyphicon-chevron-down').addClass('glyphicon-chevron-up');
	}
});

function generate_card(id) {
    var current_url = getBaseURL();
    $.getJSON(current_url + 'api/metabolites/?where={"_id" :"'+id+'"}', function(json){
        var metabolite = json["_items"][0];
        clear_card();

        $("#card_metabolite_name").html(metabolite["name"]);
        $("#molecular_formula").html(metabolite["molecular_formula"].replace(/([0-9]+)/g, '<sub>$1</sub>'));

        for (result in metabolite["adducts"]) {
            if (result == "neutral") {
                $("#neutral_mass").html(metabolite["adducts"]["neutral"][0]["accurate_mass"]);
            }

            else {
                for (adduct_idx in metabolite["adducts"][result]) {
                    var adduct = metabolite["adducts"][result][adduct_idx];
                    $("#" + result + "_adduct").append("<li class='list-group-item'><b>" + adduct["type"] + ":</b> " + adduct["accurate_mass"].toFixed(4)+"</li>");
                }
            }
        }

        if (metabolite["origins"] == null) {
            $("#origins").append("<i class='text-primary'>None</i>");
        }
        else {
            for (indx in metabolite["origins"]) {
                $("#origins").append(metabolite["origins"][indx] + "; ");
            }
        }

        if (metabolite["biofluid_locations"] == null) {
            $("#biofluids").append("<i class='text-primary'>None</i>");
        }
        else {
            for (indx in metabolite["biofluid_locations"]) {
                $("#biofluids").append(metabolite["biofluid_locations"][indx] + "; ");
            }
        }

        if (metabolite["tissue_locations"] == null) {
            $("#tissues").append("<i class='text-primary'>None</i>");
        }
        else {
            for (indx in metabolite["tissue_locations"]) {
                $("#tissues").append(metabolite["tissue_locations"][indx] + "; ");
            }
        }

        $("#viewfromcard").attr("href", current_url+"view/"+id);

    });
    $("#metabolite_card").modal("toggle");
}