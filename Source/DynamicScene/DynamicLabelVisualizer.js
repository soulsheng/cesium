/*global define*/
define([
        '../Core/Color',
        '../Core/Cartesian2',
        '../Core/Cartesian3',
        '../Scene/LabelCollection',
        '../Scene/LabelStyle',
        '../Scene/HorizontalOrigin',
        '../Scene/VerticalOrigin'
    ], function(
        Color,
        Cartesian2,
        Cartesian3,
        LabelCollection,
        LabelStyle,
        HorizontalOrigin,
        VerticalOrigin) {
    "use strict";

    function DynamicLabelVisualizer(scene, dynamicObjectCollection) {
        this._scene = scene;
        this._unusedIndexes = [];
        this._dynamicObjectCollection = undefined;

        var labelCollection = this._labelCollection = new LabelCollection();
        scene.getPrimitives().add(labelCollection);
        this.setDynamicObjectCollection(dynamicObjectCollection);
    }

    DynamicLabelVisualizer.prototype.getScene = function() {
        return this._scene;
    };

    DynamicLabelVisualizer.prototype.getDynamicObjectCollection = function() {
        return this._dynamicObjectCollection;
    };

    DynamicLabelVisualizer.prototype.setDynamicObjectCollection = function(dynamicObjectCollection) {
        var oldCollection = this._dynamicObjectCollection;
        if (oldCollection !== dynamicObjectCollection) {
            if (typeof oldCollection !== 'undefined') {
                oldCollection.objectsRemoved.removeEventListener(DynamicLabelVisualizer.prototype._onObjectsRemoved);
                this.removeAll();
            }
            this._dynamicObjectCollection = dynamicObjectCollection;
            if (typeof dynamicObjectCollection !== 'undefined') {
                dynamicObjectCollection.objectsRemoved.addEventListener(DynamicLabelVisualizer.prototype._onObjectsRemoved, this);
            }
        }
    };

    DynamicLabelVisualizer.prototype.update = function(time) {
        var dynamicObjects = this._dynamicObjectCollection.getObjects();
        for ( var i = 0, len = dynamicObjects.length; i < len; i++) {
            this.updateObject(time, dynamicObjects[i]);
        }
    };

    var show;
    var position;
    var fillColor;
    var outlineColor;
    var style;
    var text;
    var font;
    var eyeOffset;
    var pixelOffset;
    var scale;
    var verticalOrigin;
    var horizontalOrigin;
    DynamicLabelVisualizer.prototype.updateObject = function(time, dynamicObject) {
        var dynamicLabel = dynamicObject.label;
        if (typeof dynamicLabel === 'undefined') {
            return;
        }

        var textProperty = dynamicLabel.text;
        if (typeof textProperty === 'undefined') {
            return;
        }

        var positionProperty = dynamicObject.position;
        if (typeof positionProperty === 'undefined') {
            return;
        }

        var label;
        var objectId = dynamicObject.id;
        var showProperty = dynamicLabel.show;
        var labelVisualizerIndex = dynamicObject.labelVisualizerIndex;
        show = dynamicObject.isAvailable(time) && (typeof showProperty === 'undefined' || showProperty.getValue(time, show));

        if (!show) {
            //don't bother creating or updating anything else
            if (typeof labelVisualizerIndex !== 'undefined') {
                label = this._labelCollection.get(labelVisualizerIndex);
                label.setShow(false);
                this._unusedIndexes.push(labelVisualizerIndex);
                dynamicObject.labelVisualizerIndex = undefined;
            }
            return;
        }

        if (typeof labelVisualizerIndex === 'undefined') {
            var unusedIndexes = this._unusedIndexes;
            var length = unusedIndexes.length;
            if (length > 0) {
                labelVisualizerIndex = unusedIndexes.pop();
                label = this._labelCollection.get(labelVisualizerIndex);
            } else {
                labelVisualizerIndex = this._labelCollection.getLength();
                label = this._labelCollection.add();
            }
            dynamicObject.labelVisualizerIndex = labelVisualizerIndex;
            label.id = objectId;

            // CZML_TODO Determine official defaults
            label.setText('');
            label.setScale(1.0);
            label.setFont('30px sans-serif');
            label.setFillColor(Color.WHITE);
            label.setOutlineColor(Color.BLACK);
            label.setStyle(LabelStyle.FILL);
            label.setPixelOffset(Cartesian2.ZERO);
            label.setEyeOffset(Cartesian3.ZERO);
            label.setHorizontalOrigin(HorizontalOrigin.CENTER);
            label.setVerticalOrigin(VerticalOrigin.CENTER);
        } else {
            label = this._labelCollection.get(labelVisualizerIndex);
        }

        label.setShow(show);

        text = textProperty.getValue(time, text);
        if (typeof text !== 'undefined') {
            label.setText(text);
        }

        position = positionProperty.getValueCartesian(time);//, position);
        if (typeof position !== 'undefined') {
            label.setPosition(position);
        }

        var property = dynamicLabel.scale;
        if (typeof property !== 'undefined') {
            scale = property.getValue(time, scale);
            if (typeof scale !== 'undefined') {
                label.setScale(scale);
            }
        }

        property = dynamicLabel.font;
        if (typeof property !== 'undefined') {
            font = property.getValue(time, font);
            if (typeof font !== 'undefined') {
                label.setFont(font);
            }
        }

        property = dynamicLabel.fillColor;
        if (typeof property !== 'undefined') {
            fillColor = property.getValue(time, fillColor);
            if (typeof fillColor !== 'undefined') {
                label.setFillColor(fillColor);
            }
        }

        property = dynamicLabel.outlineColor;
        if (typeof property !== 'undefined') {
            outlineColor = property.getValue(time, outlineColor);
            if (typeof outlineColor !== 'undefined') {
                label.setOutlineColor(outlineColor);
            }
        }

        property = dynamicLabel.style;
        if (typeof property !== 'undefined') {
            style = property.getValue(time, style);
            if (typeof style !== 'undefined') {
                label.setStyle(style);
            }
        }

        property = dynamicLabel.pixelOffset;
        if (typeof property !== 'undefined') {
            pixelOffset = property.getValue(time, pixelOffset);
            if (typeof pixelOffset !== 'undefined') {
                label.setPixelOffset(pixelOffset);
            }
        }

        property = dynamicLabel.eyeOffset;
        if (typeof property !== 'undefined') {
            eyeOffset = property.getValue(time, eyeOffset);
            if (typeof eyeOffset !== 'undefined') {
                label.setEyeOffset(eyeOffset);
            }
        }

        property = dynamicLabel.horizontalOrigin;
        if (typeof property !== 'undefined') {
            horizontalOrigin = property.getValue(time, horizontalOrigin);
            if (typeof horizontalOrigin !== 'undefined') {
                label.setHorizontalOrigin(horizontalOrigin);
            }
        }

        property = dynamicLabel.verticalOrigin;
        if (typeof property !== 'undefined') {
            verticalOrigin = property.getValue(time, verticalOrigin);
            if (typeof verticalOrigin !== 'undefined') {
                label.setVerticalOrigin(verticalOrigin);
            }
        }
    };

    DynamicLabelVisualizer.prototype.removeAll = function() {
        this._unusedIndexes = [];
        this._labelCollection.removeAll();
        var dynamicObjects = this._dynamicObjectCollection.getObjects();
        for ( var i = dynamicObjects.length - 1; i > -1; i--) {
            dynamicObjects[i].labelVisualizerIndex = undefined;
        }
    };

    DynamicLabelVisualizer.prototype._onObjectsRemoved = function(dynamicObjectCollection, dynamicObjects) {
        var thisLabelCollection = this._labelCollection;
        var thisUnusedIndexes = this._unusedIndexes;
        for ( var i = dynamicObjects.length - 1; i > -1; i--) {
            var dynamicObject = dynamicObjects[i];
            var labelVisualizerIndex = dynamicObject.labelVisualizerIndex;
            if (typeof labelVisualizerIndex !== 'undefined') {
                var label = thisLabelCollection.get(labelVisualizerIndex);
                label.setShow(false);
                thisUnusedIndexes.push(labelVisualizerIndex);
                dynamicObject.labelVisualizerIndex = undefined;
            }
        }
    };

    return DynamicLabelVisualizer;
});