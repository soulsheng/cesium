/*global define*/
define([
        '../Shaders/Materials/DistanceIntervalMaterial'
    ], function(
        ShadersDistanceIntervalMaterial) {
    "use strict";

    /**
     * DOC_TBA
     *
     * @alias DistanceIntervalMaterial
     * @constructor
     */
    var DistanceIntervalMaterial = function(template) {
        var t = template || {};

        /**
         * DOC_TBA
         */
        this.intervals = t.intervals || [];

        // TODO: Expose get/set - can change distance/color, but not number of intervals
        var distances = [];
        var colors = [];

        for ( var i = 0; i < this.intervals.length; ++i) {
            distances.push(this.intervals[i].distance);
            colors.push(this.intervals[i].color);
        }

        this._uniforms = {
            u_distances : function() {
                return distances;
            },
            u_colors : function() {
                return colors;
            }
        };

        this.shaderSource = '#define NUMBER_OF_DISTANCES ' + this.intervals.length.toString() + '\n' +
                            '#line 0\n' +
                            ShadersDistanceIntervalMaterial;
    };

    return DistanceIntervalMaterial;
});