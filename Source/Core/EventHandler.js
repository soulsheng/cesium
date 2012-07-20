/*global define*/
define([
        './DeveloperError',
        './destroyObject',
        './Cartesian2',
        './JulianDate',
        './MouseEventType',
        './EventModifier'
    ], function(
        DeveloperError,
        destroyObject,
        Cartesian2,
        JulianDate,
        MouseEventType,
        EventModifier) {
    "use strict";

    /**
     * Handles user input events. Custom functions can be added to be executed on
     * when the user enters input.
     *
     * @alias EventHandler
     *
     * @param {DOC_TBA} element The element to add events to. Defaults to document.
     * @constructor
     */
    var EventHandler = function(element) {
        this._mouseEvents = {};
        for ( var button in MouseEventType) {
            if (MouseEventType.hasOwnProperty(button)) {
                this._mouseEvents[button] = 0;
            }
        }

        this._modifiedMouseEvents = {};
        for ( var modifier in EventModifier) {
            if (EventModifier.hasOwnProperty(modifier)) {
                this._modifiedMouseEvents[modifier] = {};
                for (button in MouseEventType) {
                    if (MouseEventType.hasOwnProperty(button)) {
                        this._modifiedMouseEvents[modifier][button] = 0;
                    }
                }
            }
        }

        this._leftMouseButtonDown = false;
        this._middleMouseButtonDown = false;
        this._rightMouseButtonDown = false;
        this._seenAnyTouchEvents = false;
        this._lastMouseX = 0;
        this._lastMouseY = 0;
        this._totalPixels = 0;

        // TODO: Revisit when doing mobile development. May need to be configurable
        // or determined based on the platform?
        this._clickPixelTolerance = 5;

        this._element = element || document;

        this._register();
    };

    EventHandler.prototype._getPosition = function(event) {
        if (this._element === document) {
            return {
                x : event.clientX,
                y : event.clientY
            };
        }

        var rect = this._element.getBoundingClientRect();
        return {
            x : event.clientX - rect.left,
            y : event.clientY - rect.top
        };
    };

    /**
     * Set a function to be executed on a mouse event.
     *
     * @memberof EventHandler
     *
     * @param {Function} action Function to be executed when the mouse event occurs.
     * @param {Enumeration} type The MouseEventType of mouse event.
     * @param {Enumeration} modifier A EventModifier key that is held when a <code>type</code>
     * event occurs.
     *
     * @exception {DeveloperError} action is required.
     * @exception {DeveloperError} type is required.
     *
     * @see EventHandler#getMouseAction
     * @see EventHandler#removeMouseAction
     */
    EventHandler.prototype.setMouseAction = function(action, type, modifier) {
        if (!action) {
            throw new DeveloperError('action is required.');
        }

        if (!type) {
            throw new DeveloperError('type is required.');
        }

        var mouseEvents;
        if (modifier && modifier.name) {
            mouseEvents = this._modifiedMouseEvents[modifier.name];
        } else {
            mouseEvents = this._mouseEvents;
        }

        if (type && type.name && mouseEvents) {
            mouseEvents[type.name] = action;
        }
    };

    /**
     * Returns the function to be executed on a mouse event.
     *
     * @memberof EventHandler
     *
     * @param {Enumeration} type The MouseEventType of mouse event.
     * @param {Enumeration} modifier A EventModifier key that is held when a <code>type</code>
     * event occurs.
     *
     * @exception {DeveloperError} type is required.
     *
     * @see EventHandler#setMouseAction
     * @see EventHandler#removeMouseAction
     */
    EventHandler.prototype.getMouseAction = function(type, modifier) {
        if (!type) {
            throw new DeveloperError('type is required.');
        }

        var mouseEvents;
        if (modifier && modifier.name) {
            mouseEvents = this._modifiedMouseEvents[modifier.name];
        } else {
            mouseEvents = this._mouseEvents;
        }

        if (type && type.name && mouseEvents) {
            return mouseEvents[type.name];
        }

        return undefined;
    };

    /**
     * Removes the function to be executed on a mouse event.
     *
     * @memberof EventHandler
     *
     * @param {Enumeration} type The MouseEventType of mouse event.
     * @param {Enumeration} modifier A EventModifier key that is held when a <code>type</code>
     * event occurs.
     *
     * @exception {DeveloperError} type is required.
     *
     * @see EventHandler#getMouseAction
     * @see EventHandler#setMouseAction
     */
    EventHandler.prototype.removeMouseAction = function(type, modifier) {
        if (!type) {
            throw new DeveloperError('type is required.');
        }

        var mouseEvents;
        if (modifier && modifier.name) {
            mouseEvents = this._modifiedMouseEvents[modifier.name];
        } else {
            mouseEvents = this._mouseEvents;
        }

        if (type && type.name && mouseEvents && mouseEvents[type.name]) {
            delete mouseEvents[type.name];
        }
    };

    EventHandler.prototype._getModifier = function(event) {
        if (event.shiftKey) {
            return EventModifier.SHIFT;
        } else if (event.ctrlKey) {
            return EventModifier.CTRL;
        } else if (event.altKey) {
            return EventModifier.ALT;
        }

        return undefined;
    };

    EventHandler.prototype._handleMouseDown = function(event) {
        var pos = this._getPosition(event);
        this._lastMouseX = pos.x;
        this._lastMouseY = pos.y;
        this._totalPixels = 0;
        if (this._seenAnyTouchEvents) {
            return;
        }

        var modifier = this._getModifier(event);
        var action;

        // IE_TODO:  On some versions of IE, the left-button is 1, and the right-button is 4.
        // See: http://www.unixpapa.com/js/mouse.html
        // This is not the case in Chrome Frame, so we are OK for now, but are there
        // constants somewhere?
        if (event.button === 0) {
            this._leftMouseButtonDown = true;
            action = this.getMouseAction(MouseEventType.LEFT_DOWN, modifier);
        } else if (event.button === 1) {
            this._middleMouseButtonDown = true;
            action = this.getMouseAction(MouseEventType.MIDDLE_DOWN, modifier);
        } else if (event.button === 2) {
            this._rightMouseButtonDown = true;
            action = this.getMouseAction(MouseEventType.RIGHT_DOWN, modifier);
        }

        if (action) {
            action({
                position : new Cartesian2(pos.x, pos.y)
            });
        }
        event.preventDefault();
    };

    EventHandler.prototype._handleMouseUp = function(event) {
        var modifier = this._getModifier(event);
        var action, clickAction;
        if (this._seenAnyTouchEvents) {
            return;
        }

        // IE_TODO:  On some versions of IE, the left-button is 1, and the right-button is 4.
        // See: http://www.unixpapa.com/js/mouse.html
        // This is not the case in Chrome Frame, so we are OK for now, but are there
        // constants somewhere?
        if (event.button === 0) {
            this._leftMouseButtonDown = false;
            action = this.getMouseAction(MouseEventType.LEFT_UP, modifier);
            clickAction = this.getMouseAction(MouseEventType.LEFT_CLICK, modifier);
        } else if (event.button === 1) {
            this._middleMouseButtonDown = false;
            action = this.getMouseAction(MouseEventType.MIDDLE_UP, modifier);
            clickAction = this.getMouseAction(MouseEventType.MIDDLE_CLICK, modifier);
        } else if (event.button === 2) {
            this._rightMouseButtonDown = false;
            action = this.getMouseAction(MouseEventType.RIGHT_UP, modifier);
            clickAction = this.getMouseAction(MouseEventType.RIGHT_CLICK, modifier);
        }

        var pos = this._getPosition(event);

        var xDiff = this._lastMouseX - pos.x;
        var yDiff = this._lastMouseY - pos.y;
        this._totalPixels += Math.sqrt(xDiff * xDiff + yDiff * yDiff);

        if (action) {
            action({
                position : new Cartesian2(pos.x, pos.y)
            });
        }

        if (clickAction && this._totalPixels < this._clickPixelTolerance) {
            clickAction({
                position : new Cartesian2(pos.x, pos.y)
            });
        }
    };

    EventHandler.prototype._handleMouseMove = function(event) {
        var pos = this._getPosition(event);
        if (this._seenAnyTouchEvents) {
            return;
        }

        var xDiff = this._lastMouseX - pos.x;
        var yDiff = this._lastMouseY - pos.y;
        this._totalPixels += Math.sqrt(xDiff * xDiff + yDiff * yDiff);

        var movement = {
            startPosition : new Cartesian2(this._lastMouseX, this._lastMouseY),
            endPosition : new Cartesian2(pos.x, pos.y),
            motion : new Cartesian2(0.0, 0.0)
        };

        var modifier = this._getModifier(event);
        var action = this.getMouseAction(MouseEventType.MOVE, modifier);
        if (action) {
            action(movement);
        }

        this._lastMouseX = movement.endPosition.x;
        this._lastMouseY = movement.endPosition.y;

        if (this._leftMouseButtonDown || this._middleMouseButtonDown || this._rightMouseButtonDown) {
            event.preventDefault();
        }
    };

    EventHandler.prototype._handleTouchStart = function(event) {
        var pos, numberOfTouches = event.touches.length;
        this._seenAnyTouchEvents = true;

        if (numberOfTouches === 1) {
            pos = this._getPosition(event.touches[0]);
            this._lastMouseX = pos.x;
            this._lastMouseY = pos.y;
            this._totalPixels = 0;

            var modifier = this._getModifier(event);
            var action;

            this._leftMouseButtonDown = true;
            action = this.getMouseAction(MouseEventType.LEFT_DOWN, modifier);

            if (action) {
                action({
                    position : new Cartesian2(pos.x, pos.y)
                });
            }
            event.preventDefault();
        } else if (this._leftMouseButtonDown) {
            this._handleTouchEnd(event);
        }
    };

    EventHandler.prototype._handleTouchEnd = function(event) {
        var numberOfTouches = event.touches.length;
        var numberOfTargetTouches = event.targetTouches.length;
        var modifier = this._getModifier(event);
        var action, clickAction;

        if (this._leftMouseButtonDown) {
            this._leftMouseButtonDown = false;
            action = this.getMouseAction(MouseEventType.LEFT_UP, modifier);
            clickAction = this.getMouseAction(MouseEventType.LEFT_CLICK, modifier);
        }

        if (numberOfTargetTouches > 0) {
            var pos = this._getPosition(event.targetTouches[0]);

            var xDiff = this._lastMouseX - pos.x;
            var yDiff = this._lastMouseY - pos.y;
            this._totalPixels += Math.sqrt(xDiff * xDiff + yDiff * yDiff);

            if (action) {
                action({
                    position : new Cartesian2(pos.x, pos.y)
                });
            }

            if (clickAction && this._totalPixels < this._clickPixelTolerance) {
                clickAction({
                    position : new Cartesian2(pos.x, pos.y)
                });
            }
        }

        if (numberOfTouches === 1) {
            this._handleTouchStart(event);
        }
    };

    EventHandler.prototype._handleTouchMove = function(event) {
        if (this._leftMouseButtonDown && (event.touches.length === 1)) {
            var pos = this._getPosition(event.touches[0]);

            var xDiff = this._lastMouseX - pos.x;
            var yDiff = this._lastMouseY - pos.y;
            this._totalPixels += Math.sqrt(xDiff * xDiff + yDiff * yDiff);

            var movement = {
                startPosition : new Cartesian2(this._lastMouseX, this._lastMouseY),
                endPosition : new Cartesian2(pos.x, pos.y),
                motion : new Cartesian2(0.0, 0.0)
            };

            var modifier = this._getModifier(event);
            var action = this.getMouseAction(MouseEventType.MOVE, modifier);
            if (action) {
                action(movement);
            }

            this._lastMouseX = movement.endPosition.x;
            this._lastMouseY = movement.endPosition.y;

            if (this._leftMouseButtonDown || this._middleMouseButtonDown || this._rightMouseButtonDown) {
                event.preventDefault();
            }
        }
    };

    EventHandler.prototype._handleMouseWheel = function(event) {
        // Some browsers use event.detail to count the number of clicks. The sign
        // of the integer is the direction the wheel is scrolled. In that case, convert
        // to the angle it was rotated in degrees.
        var delta = event.detail ? event.detail * -120 : event.wheelDelta;

        var modifier = this._getModifier(event);
        var type = MouseEventType.WHEEL;
        var action = this.getMouseAction(type, modifier);

        if (action) {
            event.preventDefault();
            action(delta);
        }
    };

    EventHandler.prototype._handleMouseDblClick = function(event) {
        var modifier = this._getModifier(event);
        var action;
        var pos = this._getPosition(event);

        // IE_TODO:  On some versions of IE, the left-button is 1, and the right-button is 4.
        // See: http://www.unixpapa.com/js/mouse.html
        // This is not the case in Chrome Frame, so we are OK for now, but are there
        // constants somewhere?
        if (event.button === 0) {
            action = this.getMouseAction(MouseEventType.LEFT_DOUBLE_CLICK, modifier);
        } else if (event.button === 1) {
            action = this.getMouseAction(MouseEventType.MIDDLE_DOUBLE_CLICK, modifier);
        } else if (event.button === 2) {
            action = this.getMouseAction(MouseEventType.RIGHT_DOUBLE_CLICK, modifier);
        }

        if (action) {
            action({
                position : new Cartesian2(pos.x, pos.y)
            });
        }
    };

    EventHandler.prototype._register = function() {
        var that = this, useDoc = true;

        this._callbacks = [];
        if (typeof this._element.disableRootEvents !== 'undefined') {
            useDoc = false;
        }

        this._callbacks.push({
            name : 'mousedown',
            onDoc : false,
            action : function(e) {
                that._handleMouseDown(e);
            }
        });
        this._callbacks.push({
            name : 'mouseup',
            onDoc : useDoc,
            action : function(e) {
                that._handleMouseUp(e);
            }
        });
        this._callbacks.push({
            name : 'mousemove',
            onDoc : useDoc,
            action : function(e) {
                that._handleMouseMove(e);
            }
        });
        this._callbacks.push({
            name : 'dblclick',
            onDoc : false,
            action : function(e) {
                that._handleMouseDblClick(e);
            }
        });
        this._callbacks.push({
            name : 'touchstart',
            onDoc : false,
            action : function(e) {
                that._handleTouchStart(e);
            }
        });
        this._callbacks.push({
            name : 'touchend',
            onDoc : useDoc,
            action : function(e) {
                that._handleTouchEnd(e);
            }
        });
        this._callbacks.push({
            name : 'touchmove',
            onDoc : useDoc,
            action : function(e) {
                that._handleTouchMove(e);
            }
        });

        // Firefox calls the mouse wheel event 'DOMMouseScroll', all others use 'mousewheel'
        this._callbacks.push({
            name : 'mousewheel',
            onDoc : false,
            action : function(e) {
                that._handleMouseWheel(e);
            }
        });
        this._callbacks.push({
            name : 'DOMMouseScroll',
            onDoc : false,
            action : function(e) {
                that._handleMouseWheel(e);
            }
        });

        for ( var i = 0; i < this._callbacks.length; i++) {
            var cback = this._callbacks[i];
            if (cback.onDoc) {
                document.addEventListener(cback.name, cback.action, false);
            } else {
                this._element.addEventListener(cback.name, cback.action, false);
            }
        }
    };

    EventHandler.prototype._unregister = function() {
        for ( var i = 0; i < this._callbacks.length; i++) {
            var cback = this._callbacks[i];
            if (cback.onDoc) {
                document.removeEventListener(cback.name, cback.action, false);
            } else {
                this._element.removeEventListener(cback.name, cback.action, false);
            }
        }
    };

    /**
     * Returns true if this object was destroyed; otherwise, false.
     * <br /><br />
     * If this object was destroyed, it should not be used; calling any function other than
     * <code>isDestroyed</code> will result in a {@link DeveloperError} exception.
     *
     * @memberof EventHandler
     *
     * @return {Boolean} <code>true</code> if this object was destroyed; otherwise, <code>false</code>.
     *
     * @see EventHandler#destroy
     */
    EventHandler.prototype.isDestroyed = function() {
        return false;
    };

    /**
     * Removes mouse listeners held by this object.
     * <br /><br />
     * Once an object is destroyed, it should not be used; calling any function other than
     * <code>isDestroyed</code> will result in a {@link DeveloperError} exception.  Therefore,
     * assign the return value (<code>undefined</code>) to the object as done in the example.
     *
     * @memberof EventHandler
     *
     * @return {undefined}
     *
     * @exception {DeveloperError} This object was destroyed, i.e., destroy() was called.
     *
     * @see EventHandler#isDestroyed
     *
     * @example
     * handler = handler && handler.destroy();
     */
    EventHandler.prototype.destroy = function() {
        this._unregister();
        return destroyObject(this);
    };

    return EventHandler;
});
