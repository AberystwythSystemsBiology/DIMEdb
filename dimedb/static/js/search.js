$(document).ready(function () {
    $("input[type=radio][name=ionisation]").change(function() {

        var iad = {
            "negative" : [["M-H", 1], ["M+Cl", 1], ["M+Br", 0]],
            "positive" : [["M+H", 1], ["M+K", 1], ["2M+H", 0]],
            "neutral" : [["M", 1]]
        }

        $("#ionisation-adducts").empty()

        for(x in iad[this.value]) {
            var adduct = (iad[this.value][x]);
            if (adduct[1] == 1) {
                $("#ionisation-adducts").append("<option selected>"+adduct[0]+"</option>")
            }
            else {
                $("#ionisation-adducts").append("<option>"+adduct[0]+"</option>")
            }
        }

    });

    $("#search_button").click(function() {
        var ionisation = $("input[type=radio][name=ionisation]:checked").val();
        var current_url = window.location.protocol + "//" + window.location.host

        var base_api = "/api/adducts/?adducts__pi_ppm=";
        base_api += ionisation + ",";
        base_api += $("#mass").val() +",";
        base_api += $("#ppm_tolerance").val();
        // console.log(current_url+base_api)

        $("input[name='origins']:checked").each(function() {
            base_api += "&origins__contains="+this.value +",";
        });

        $("input[name='biofluids']:checked").each(function() {
            base_api += "&biofluids__contains="+this.value +",";
        });

        console.log(current_url+base_api)
    });

});

// $("#ionisation-adducts").val()
